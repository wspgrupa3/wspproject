#install.packages("gWidgetsRGtk2", dep = TRUE)
rm(list = ls());
options(show.error.locations=TRUE);
require(gWidgets)
options(guiToolkit = "RGtk2")
set_wd_click <- function(h, ...){#loading files and GUI when wd known
  user_wd = svalue(h$action);
  tryCatch(
    {
      user_wd = gsub("\\\\","/", user_wd);
      setwd(user_wd);
      source("R_WSP_gui_functions.R");#ok
      source("get_expset.R");#ok
      source("gene_diff.R");#lille error
      source("Klasteryzacja_MS.R");#ok
      source("klasteryzacja.R");#error commented out; no output
      source("klas_hier.R");#?
      addHandlerClicked(wd_select_file_button, handler = get_file, action = selfile_label);
      addHandlerClicked(process_raw_dat_button, handler = preprocess_CEL_data);
      addHandlerClicked(analysis_pick_button, handler = create_tools);
    },error = function(e){
      svalue(h$action) = e;#'Wrong directory';
    }
  )
}
#basic GUI creation
#items in gframes, frames in 3 main groups horizontally
appwin <- gwindow();
side_container = ggroup(horizontal = F, container = appwin);
wdf = gframe(text = "Working direcory", horizontal = F, container = side_container);
wdfu = ggroup(horizontal = T);
wd_paste = gedit("Paste_here", container = wdfu);
wd_select_button = gbutton("Set wd", container = wdfu);
addHandlerClicked(wd_select_button, handler=set_wd_click, action = wd_paste);
wdfd = ggroup(horizontal = T);
selfile_label = glabel(text = "NO FILE", container = wdfd);
wd_select_file_button = gbutton("Choose file", container = wdfd);
add(wdf, wdfu);
add(wdf, wdfd);
process_raw_dat_button = gbutton("Preprocess dataset", container = wdf);
asf = gframe(text = "Analysys options", container = side_container, horizontal = F);
analysis_radio = gradio(c("Genes Difference", "Clustering", "Clustering 2", "PCA"), container = asf);
analysis_pick_button = gbutton("Choose method", container = asf);
settings_container = gframe(text = "Settings", horizontal = F, container = appwin);
size(settings_container) = c(200,300);
results_container = ggroup(horizontal = F, container = appwin);
plot_frame = gframe(text = " Results", horizontal = T, container = results_container);
ggraphics(container = plot_frame);
plot_device_1 = dev.cur();
ggraphics(container = plot_frame);
plot_device_2 = dev.cur();
results_text = glabel(container = results_container);
