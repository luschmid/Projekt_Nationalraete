path_local <- "C:/Schmidlu/Dropbox/Projekt Nationalräte/04_Results/04_Political_Rents"
path_overleaf <- "C:/Schmidlu/Dropbox/Apps/Overleaf/Political_Rents"
setwd(path_local)

# A) Functions

FigureOutExt <- function(input,filename){
return(
cat(paste0("\\begin{figure}[!ht]
\\begin{center}
\\caption{Extensive margin}\\vspace{0.5em}
\\begin{tabular}{cccc}
\\includegraphics[scale=0.3]{figures/fig_",input,"_i_all_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_all_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_all_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_i_all_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_i_lrg_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_lrg_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_lrg_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_i_lrg_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_i_sml_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_sml_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_sml_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_i_sml_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_i_prs_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_prs_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_i_prs_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_i_prs_s2.pdf}\\\\
\\end{tabular}
\\end{center}
\\end{figure}"),file=filename, append=TRUE)
)}

FigureOutIntSum <- function(input,filename){
return(
cat(paste0("
\\begin{figure}[!ht]
\\begin{center}
\\caption{Intensive margin (sum)}
\\vspace{0.5em}
\\begin{tabular}{cccc}
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_sum_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_sum_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_sum_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_sum_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_sum_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_sum_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_sum_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_sum_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_sum_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_sum_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_sum_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_sum_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_sum_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_sum_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_sum_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_sum_s2.pdf}\\\\
\\end{tabular}
\\end{center}
\\end{figure}"),file=filename, append=TRUE))
}


FigureOutIntAvg <- function(input,filename){
  return(
    cat(paste0("
\\begin{figure}[!ht]
\\begin{center}
\\caption{Intensive margin (avg)}
\\vspace{0.5em}
\\begin{tabular}{cccc}
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_avg_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_avg_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_avg_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_all_avg_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_avg_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_avg_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_avg_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_lrg_avg_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_avg_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_avg_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_avg_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_sml_avg_s2.pdf}\\\\
\\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_avg_c1.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_avg_c2.pdf} & \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_avg_s1.pdf}& \\includegraphics[scale=0.3]{figures/fig_",input,"_n_prs_avg_s2.pdf}\\\\
\\end{tabular}
\\end{center}
\\end{figure}
\\clearpage"),file=filename, append=TRUE))
}


FigureOutRound2 <- function(input,filename){
  return(
    cat(paste0("
\\begin{figure}[!ht]
\\begin{center}
\\caption{}
\\vspace{0.5em}
\\begin{tabular}{cccc}
\\includegraphics[scale=0.45]{figures/fig_",input,"_i_all_c1.pdf} & \\includegraphics[scale=0.45]{figures/fig_",input,"_n_all_sum_c1.pdf}  \\\\
\\includegraphics[scale=0.45]{figures/fig_",input,"_i_lrg_c1.pdf} & \\includegraphics[scale=0.45]{figures/fig_",input,"_n_lrg_sum_c1.pdf}  \\\\
\\includegraphics[scale=0.45]{figures/fig_",input,"_i_sml_c1.pdf}  & \\includegraphics[scale=0.45]{figures/fig_",input,"_n_sml_sum_c1.pdf} \\\\
 \\includegraphics[scale=0.45]{figures/fig_",input,"_i_prs_c1.pdf} & \\includegraphics[scale=0.45]{figures/fig_",input,"_n_prs_sum_c1.pdf} \\\\
\\end{tabular}
\\end{center}
\\end{figure}
\\clearpage"),file=filename, append=TRUE))
}


FigureOutRound2_Means <- function(input,filename){
outfiles <-  paste0("\\includegraphics[scale=0.45]{figures/",input[1],"_ob.pdf} & \\includegraphics[scale=0.45]{figures/",input[1],"_hb.pdf}  \\\\")
for (i in c(2:length(input))){
outfiles <- paste0(outfiles, 
                  "\\includegraphics[scale=0.45]{figures/",input[i],"_ob.pdf} & \\includegraphics[scale=0.45]{figures/",input[i],"_hb.pdf}  \\\\")
}  
  
out <- paste0("
\\begin{figure}[!ht]
\\begin{center}
\\caption{}
\\vspace{0.5em}
\\begin{tabular}{cccc}",
outfiles,
"\\end{tabular}
\\end{center}
\\end{figure}
\\clearpage")
return(cat(out,file=filename, append=FALSE))
}


# B) Round 1: Loop over all files
# Note: Discussed by Mark, Simon, and Lukas on June 24, 2024 in Lucerne. 

es <- c("rb","cl")
bw <- c("ob","hb")
po <- c("p1","p2")
ke <- c("tri","epa","uni")
ye <- c("ay","1y")

x <- expand.grid(es,bw,po,ke,ye, stringsAsFactors = FALSE)

loop_files <- paste(x$Var1,x$Var2,x$Var3,x$Var4,x$Var5,sep="_")

#loop_files1 <- loop_files[1:(length(loop_files)/2)]
#loop_files2 <- loop_files[(length(loop_files)/2+1):length(loop_files)]

try(file.remove("out.txt"))

for (f in loop_files){
  FigureOutExt(input=f,file="out_round1.txt")
  FigureOutIntSum(input=f,file="out_round1.txt")
  FigureOutIntAvg(input=f,file="out_round1.txt")
}


# C) Round 2: Loop over selected files

es <- c("rb")
bw <- c("ob","hb")
po <- c("p1","p2")
ke <- c("tri")
ye <- c("ay","1y")

x <- expand.grid(es,bw,po,ke,ye, stringsAsFactors = FALSE)

loop_files <- paste(x$Var1,x$Var2,x$Var3,x$Var4,x$Var5,sep="_")

try(file.remove("out_round2.txt"))

for (f in loop_files){
  FigureOutRound2(input=f,file="out_round2.txt")
}


# D) Round 2: Mean comparison of samples


var <- c("year","birthyear","canton_lrg","leftist","centrist","rightist")
loop_files <- paste("fig_mean_i_lrg_c1",var,sep="_")

FigureOutRound2_Means(input=loop_files[1:3],filename="out_round2_mean_comparison1.txt")
FigureOutRound2_Means(input=loop_files[4:6],filename="out_round2_mean_comparison2.txt")




