#install.packages("gWidgetsRGtk2", dep = TRUE)
require(gWidgets)
options(guiToolkit = "RGtk2")
set_wd_click <- function(h, ...){#loading files and GUI when wd known
  user_wd = svalue(h$action);
  tryCatch(
    {
      user_wd = gsub("\\\\","/", user_wd);
      setwd(user_wd);
      source("R_WSP_gui_functions.R");
      source("get_expset.R");
      source("klasteryzacja.R");
      source("klas_hier.R");
      source("pca.R");
      gene_diff_group_input = gframe(text = "Genes difference", horizontal = F);
      glabel(text = "FC", container = gene_diff_group_input);
      pval_slider = gslider(from = 0.01, to = 0.1, by = 0.01, value = 0.05, container = gene_diff_group_input);
      glabel(text = "p-val threshold", container = gene_diff_group_input);
      FC_slider = gslider(from = -7, to = 7, by = 0.5, value = 0, container = gene_diff_group_input);
      add(central_container, gene_diff_group_input);
      gene_diff_group_output = ggroup(horizontal = F);
      #pca
      pca_group_input = gframe(text = "PCA", horizontal = F);
      glabel("Number of clusters", container = pca_group_input);
      clnum_slider = gslider(from = 2, to = 10, by = 1, value = 3, container = pca_group_input);
      add(central_container, pca_group_input);
      pca_group_output = gframe("PCA plots", horizontal = T);
      #pca_win1 = ggraphics(container = pca_group_output);
      #pca_win2 = ggraphics(container = pca_group_output);
      #pca_win3 = ggraphics(container = pca_group_output);
      add(results_container, pca_group_output);
      addHandlerClicked(wd_select_file_button, handler = get_file, action = selfile_label);
      addHandlerClicked(process_raw_dat_button, handler = process_file);
    },error = function(e){
      svalue(h$action) = 'Wrong directory';
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
asf = gframe(text = "Analysys options", container = side_container);
analysis_radio = gradio(c("F1", "F2"), container = asf);
central_container = gframe(text = "Settings", horizontal = F, container = appwin);
results_container = ggroup(horizontal = F, container = appwin);
plot_device = ggraphics(container = results_container);