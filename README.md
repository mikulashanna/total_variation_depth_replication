# Partial Replication of Huang Huang, Ying Sun: A Decomposition of Total Variation Depth for Understanding Functional Outliers (2019)

## Contents
* [observations_and_questions.Rmd](observations_and_questions.Rmd) containts my observations and questions about the paper.
* Within [src](src),
  + [data_generation_functions.R](src/data_generation_functions.R) containts functions that generate data,
  + [analysis_functions.R](src/analysis_functions.R) containts functions that analyze data, and
  + [visualization_functions.R](src/visualization_functions.R) containts functions that generate figures or tables.
* [data](data) contains data I used to replicate figures and tables.
* [generate_output.Rmd](generate_output.Rmd) replicates figures and tables, updates the contents of [output](output), and records the progress for code chunks `table 1 chunk 1`, `table 1 chunk 2`, and `table 1 test` at [progress](progress).
* Within [docs](docs),
  + [observations_and_questions.html](docs/observations_and_questions.html) is the `.html` file that `rmarkdown::render(input="observations_and_questions.Rmd",output_file="docs/observations_and_questions.html")` created, and
  + [generate_output.html](docs/generate_output.html) is the `.html` file that `rmarkdown::render(input="generate_output.Rmd",output_file="docs/generate_output.html")` created.
* [output](output) contains the replicated figures and tables.
* [progress](progress) contains the progress bars for code chunks `table 1 chunk 1`, `table 1 chunk 2`, and `table 1 test` of [generate_output.Rmd](generate_output.Rmd).

The rendered version of [docs/observations_and_questions.html](docs/observations_and_questions.html) is at <https://mikulashanna.github.io/total_variation_depth_replication/observations_and_questions.html>, and the rendered version of [docs/generate_output.html](docs/generate_output.html) is at <https://mikulashanna.github.io/total_variation_depth_replication/generate_output.html>.

## Reproducibility

Line 37 of [generate_output.Rmd](generate_output.Rmd) sets the seed.

## References
* original paper:<br>
  + Huang, H., &amp; Sun, Y. (2019). A decomposition of total variation depth for understanding functional outliers. *Technometrics, 61*(4), 445–458. <https://doi.org/10.1080/00401706.2019.1574241>
* original code:<br>
  + Huang, H., &amp; Sun, Y. (2018, July 28). *TVD.R.* GitHub. <https://github.com/hhuang90/TVD/blob/master/R/TVD.R>
* source of [data/figure_5_data.txt](data/figure_5_data.txt):<br>
  + Climate Prediction Center Internet Team. (2021, January 21). *Monthly Atmospheric &amp; Sea Surface Temperature Indices.* Climate Prediction Center. <https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices>
* papers used for creating functions in [src](src):<br>
  + López-Pintado, S., &amp; Romo, J. (2009). On the concept of depth for functional data. *Journal of the American Statistical Association, 104*(486), 718–734. <https://doi.org/10.1198/jasa.2009.0108><br>
  + Sun, Y., &amp; Genton, M. G. (2011). Functional boxplots. *Journal of Computational and Graphical Statistics, 20*(2), 316–334. https://doi.org/10.1198/jcgs.2011.09224<br>
  + Narisetty, N. N., &amp; Nair, V. N. (2016). Extremal depth for functional data and applications. *Journal of the American Statistical Association, 111*(516), 1705–1714. https://doi.org/10.1080/01621459.2015.1110033
