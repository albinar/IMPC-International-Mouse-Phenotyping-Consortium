This assumes that you have your results from the pipeline.
You need to set your input directory which is you results. this folder contain the follwing folder:
+flowcut
++flowcut/plots
+gating (there can be multiple folders here, just make sure your automatedplots have the suffix of "_Automated.png"
+stats
++stats/counts (there can be multiple folders here, just make sure your csv files have the suffix of "_counts.csv"
+workspace



In ui.R, you need to change the following:
setwd() path
titles
any of the buttons, if you want to label them differently

In server.R, you need to make sure your csv files are stored in a way explained in shinyHelper.R, otherwise you need to change that function and corresponding outputs accordingly

