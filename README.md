**AIConsult.LC**

AIConsult.LC is an R package designed to estimate individualized lung cancer risk predictions using globally validated risk models. It provides tools for batch risk computation and interactive visualization of results, tailored for researchers, clinicians and even people in general.

**Key Features:**
1. Batch Prediction: Compute lung cancer risks for multiple individuals via lcrisks().
2. Individual Visualization: Generate tailored plots for single subjects using individualPred_plot().

**Installation**

1. Install from Local File：Download the AIConsult.LC_0.1.0.tar.gz file from the Releases page.
2. Install in R/RStudio: Navigate to: Tools → Install Packages → Select Install from Package Archive File → Browse to the downloaded .tar.gz file and install.

**Reference models**
1. Model: Bach
   
   Reference：Bach PB, Kattan MW, Thornquist MD, et al. Variations in lung cancer risk among smokers. J Natl Cancer Inst. 2003;95(6):470-478. doi:10.1093/jnci/95.6.470
2. Model: Spitz

   Reference：Spitz MR, Hong WK, Amos CI, et al. A risk model for prediction of lung cancer. J Natl Cancer Inst. 2007;99(9):715-726. doi:10.1093/jnci/djk153
3. Model: Hoggart
   
   Reference：Hoggart C, Brennan P, Tjonneland A, et al. A risk model for lung cancer incidence. Cancer Prev Res (Phila). 2012;5(6):834-846. doi:10.1158/1940-6207.CAPR-11-0237
4. Model: PLCOm2012
   
   Reference：Tammemägi MC, Katki HA, Hocking WG, et al. Selection criteria for lung-cancer screening [published correction appears in N Engl J Med. 2013 Jul 25;369(4):394]. N Engl J Med. 2013;368(8):728-736. doi:10.1056/NEJMoa1211776
5. Model: Korean Men
   
   Reference：Park S, Nam BH, Yang HR, et al. Individualized risk prediction model for lung cancer in Korean men. PLoS One. 2013;8(2):e54823. doi:10.1371/journal.pone.0054823
6. Model: PLCOall2014
    
   Reference：Tammemägi MC, Church TR, Hocking WG, et al. Evaluation of the lung cancer risks at which to screen ever- and never-smokers: screening rules applied to the PLCO and NLST cohorts [published correction appears in PLoS Med. 2015 Jan 28;12(1):e1001787. doi: 10.1371/journal.pmed.1001787]. PLoS Med. 2014;11(12):e1001764. Published 2014 Dec 2. doi:10.1371/journal.pmed.1001764
7. Model: Pittsburgh Predictor
    
   Reference：Wilson DO, Weissfeld J. A simple model for predicting lung cancer occurrence in a lung cancer screening program: The Pittsburgh Predictor. Lung Cancer. 2015;89(1):31-37. doi:10.1016/j.lungcan.2015.03.021
8. Model: LLPi

   Reference：Marcus MW, Chen Y, Raji OY, Duffy SW, Field JK. LLPi: Liverpool Lung Project Risk Prediction Model for Lung Cancer Incidence. Cancer Prev Res (Phila). 2015;8(6):570-575. doi:10.1158/1940-6207.CAPR-14-0438
9. Model: LCRAT

   Reference：Katki HA, Kovalchik SA, Berg CD, Cheung LC, Chaturvedi AK. Development and Validation of Risk Models to Select Ever-Smokers for CT Lung Cancer Screening. JAMA. 2016;315(21):2300-2311. doi:10.1001/jama.2016.6255
10. Model: HUNT

    Reference：Markaki M, Tsamardinos I, Langhammer A, Lagani V, Hveem K, Røe OD. A Validated Clinical Risk Prediction Model for Lung Cancer in Smokers of All Ages and Exposure Types: A HUNT Study [published correction appears in EBioMedicine. 2022 Aug;82:104187. doi: 10.1016/j.ebiom.2022.104187]. EBioMedicine. 2018;31:36-46. doi:10.1016/j.ebiom.2018.03.027
11. Model: LLPv3
    
    Reference：Field JK, Vulkan D, Davies MPA, Duffy SW, Gabe R. Liverpool Lung Project lung cancer risk stratification model: calibration and prospective validation. Thorax. 2021;76(2):161-168. doi:10.1136/thoraxjnl-2020-215158
12. Model: LCRS
    
    Reference：Ma Z, Lv J, Zhu M, et al. Lung cancer risk score for ever and never smokers in China. Cancer Commun (Lond). 2023;43(8):877-895. doi:10.1002/cac2.12463
13. Model: OWL
    
    Reference：Pan Z, Zhang R, Shen S, et al. OWL: an optimized and independently validated machine learning prediction model for lung cancer screening based on the UK Biobank, PLCO, and NLST populations. EBioMedicine. 2023;88:104443. doi:10.1016/j.ebiom.2023.104443

**How to cite**

If you use this R package in your work, please cite the following foundational paper:

Ziqing Ye, Yexiang Sun, Yueqi Yin, et al. Assessment and Recalibration of Seventeen Lung Cancer Risk Prediction Models in approximately One Million Chinese Population Utilizing Healthcare Big Data: a Retrospective Cohort Analysis. The Lancet Regional Health – Western Pacific 2025.
