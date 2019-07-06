for (px in 1:sys.nframe()) {
  if (exists("itemModelExplorer.default", envir=parent.frame(px))) {
    itemModelExplorer.default <-
      get("itemModelExplorer.default", envir=parent.frame(px))
    break
  }
}

if (!exists("itemModelExplorer.default")) {
  itemModelExplorer.default <-
    list(scaleValue = 1,
         thresholds = c(.8, 1.6),
         discrimination = 1)
}
