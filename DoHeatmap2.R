##Modified Doheatmap from Seurat for aestectics - specific to singlecell_dmel_hemocytes_jul2020.R
##lines 165-167 were modified

DoHeatmap2 <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      cols <- group.colors[1:length(x = levels(x = group.use))] %||% default.colors
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
          x = cols,
          start = 1,
          stop = 7
        )))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
            x = cols,
            start = 1,
            stop = 7
          )))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      # scale the height of the bar
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      plot <- plot +
        annotation_raster(
          raster = t(x = cols[group.use2]),
          xmin = -Inf,
          xmax = Inf,
          ymin = y.pos,
          ymax = y.max
        ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off') +
        scale_color_manual(values = cols)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.5,
          angle = angle,
          hjust = hjust,
          size = size,
          vjust = 0.1,
          color= c('#9970ab',"#F28848","#FCBC2A", "#F0F921"),
          fontface = "bold"
        )
        plot <- suppressMessages(plot + coord_cartesian(
          ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
          clip = 'off')
        )
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}


`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}

CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}


SingleRasterMap <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL
) {
  data <- MinMax(data = data, min = disp.min, max = disp.max)
  data <- Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  plot <- ggplot(data = data) +
    my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors, na.value = "white") +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    WhiteBackground() + NoAxes(keep.text = TRUE)
  if (!is.null(x = group.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
      alpha = 0
    ) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  }
  return(plot)
}

Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}
