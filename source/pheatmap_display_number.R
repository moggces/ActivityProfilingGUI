pheatmap_new_label <- function (mat, mat_num, color = colorRampPalette(rev(c("#D73027", "#FC8D59", 
    "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100), 
    kmeans_k = NA, breaks = NA, border_color = "grey60", cellwidth = NA, 
    cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
    clustering_method = "complete", treeheight_row = ifelse(cluster_rows, 
        50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), 
    legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation = NA, 
    annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, 
    show_rownames = T, show_colnames = T, main = NA, fontsize = 10, 
    fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, 
    number_format = "%.2f", fontsize_number = 0.8 * fontsize, 
    filename = NA, width = NA, height = NA, ...) 
{
    mat = as.matrix(mat)
    mat = scale_mat(mat, scale)
    if (!is.na(kmeans_k)) {
        km = kmeans(mat, kmeans_k, iter.max = 100)
        mat = km$centers
        t = table(km$cluster)
        rownames(mat) = sprintf("cl%s_size_%d", names(t), t)
    }
    else {
        km = NA
    }
    if (cluster_rows) {
        tree_row = cluster_mat(mat, distance = clustering_distance_rows, 
            method = clustering_method)
        mat = mat[tree_row$order, ]
		mat_num = mat_num[tree_row$order, ]
    }
    else {
        tree_row = NA
        treeheight_row = 0
    }
    if (cluster_cols) {
        tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, 
            method = clustering_method)
        mat = mat[, tree_col$order]
		mat_num = mat_num[, tree_col$order]
    }
    else {
        tree_col = NA
        treeheight_col = 0
    }
    if (display_numbers) {
        #fmat = matrix(sprintf(number_format, mat_num), nrow = nrow(mat_num), 
        #    ncol = ncol(mat_num))
		fmat = mat_num
        attr(fmat, "draw") = TRUE
    }
    else {
        fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
        attr(fmat, "draw") = FALSE
    }
    if (!is.na(legend_breaks[1]) & !is.na(legend_labels[1])) {
        if (length(legend_breaks) != length(legend_labels)) {
            stop("Lengths of legend_breaks and legend_labels must be the same")
        }
    }
    if (is.na(breaks[1])) {
        breaks = generate_breaks(as.vector(mat), length(color))
    }
    if (legend & is.na(legend_breaks[1])) {
        legend = grid.pretty(range(as.vector(breaks)))
        names(legend) = legend
    }
    else if (legend & !is.na(legend_breaks[1])) {
        legend = legend_breaks[legend_breaks >= min(breaks) & 
            legend_breaks <= max(breaks)]
        if (!is.na(legend_labels[1])) {
            legend_labels = legend_labels[legend_breaks >= min(breaks) & 
                legend_breaks <= max(breaks)]
            names(legend) = legend_labels
        }
        else {
            names(legend) = legend
        }
    }
    else {
        legend = NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)
    if (!is.na(annotation[[1]][1])) {
        annotation = annotation[colnames(mat), , drop = F]
        annotation_colors = generate_annotation_colours(annotation, 
            annotation_colors, drop = drop_levels)
    }
    if (!show_rownames) {
        rownames(mat) = NULL
    }
    if (!show_colnames) {
        colnames(mat) = NULL
    }
    heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, 
        cellheight = cellheight, treeheight_col = treeheight_col, 
        treeheight_row = treeheight_row, tree_col = tree_col, 
        tree_row = tree_row, filename = filename, width = width, 
        height = height, breaks = breaks, color = color, legend = legend, 
        annotation = annotation, annotation_colors = annotation_colors, 
        annotation_legend = annotation_legend, main = main, fontsize = fontsize, 
        fontsize_row = fontsize_row, fontsize_col = fontsize_col, 
        fmat = fmat, fontsize_number = fontsize_number, ...)
    invisible(list(tree_row = tree_row, tree_col = tree_col, 
        kmeans = km))
}