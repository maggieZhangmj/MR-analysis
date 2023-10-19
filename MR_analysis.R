################R包自带数据检验成功#########

#准备工作

install.packages("remotes") remotes::install_github("MRCIEU/TwoSampleMR") remotes::install_github("MRCIEU/TwoSampleMR@0.4.26") library(TwoSampleMR)

#第二步：获取工具变量-暴露因素的GWAS数据

exposure_dat <-extract_instruments("ieu-a-2")

#第三步：获取暴露因素-结局变量的GWAS数据

outcome_dat <-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-7")

#第四步：对数据进行预处理，使其效应等位与效应量保持统一

dat <- harmonise_data( exposure_dat = exposure_dat, outcome_dat = outcome_dat )

#第五步：MR分析

res <- mr(dat) res

#第六步：敏感性分析（有三部分组成|结果出来需要检验是否可靠）

#第一部分：异质性检验 mr_heterogeneity

mr_heterogeneity(dat)

#第二部分：水平多效性检验Horizontal pleiotropy

mr_pleiotropy_test(dat)

#第三部分：Leave-one-out analysis

res_loo <- mr_leaveoneout(dat) mr_leaveoneout_plot(res_loo) p1 <- mr_scatter_plot(res, dat) p1

#第七步：可视化结果（散点图、森林图、漏斗图）

#森林图可视化

res_single <- mr_singlesnp(dat) res_single mr_forest_plot(res_single)

#绘制漏斗图

mr_funnel_plot(res_single)

#如果要去除点，我们需要知道，可视化的本质是数据，所以只需要查看可视化数据即可

p1 <- mr_funnel_plot(res_single) p1 ###############################################
