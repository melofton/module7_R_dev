# Macrosystems EDDIE Module 7: Using Data to Improve Ecological Forecasts

## Summary
**Ecological forecasting** is a tool that can be used for understanding and predicting changes in populations, communities, and ecosystems. Ecological forecasting is an emerging approach which provides an estimate of the future state of an ecological system with uncertainty, allowing society to prepare for changes in important ecosystem services. 

**How can we use data to improve ecological forecasts?**

To be useful for management, ecological forecasts need to be both accurate enough for managers to be able to rely on them for decision-making and include a representation of forecast uncertainty, so managers can properly interpret the probability of future events. To improve forecast accuracy, we can update forecasts with observational data once they become available, a process known as **data assimilation**. Recent improvements in environmental sensor technology and an increase in the number of sensors deployed in ecosystems have increased the availability of data for assimilation to develop and improve forecasts for natural resource management. 

In this module, you will explore how assimilating data with different amounts of observation uncertainty and at different temporal frequencies affects forecasts of water quality at a lake site of your choice. 
  
## Learning Outcomes
1. Define data assimilation. (Activity A)    
2. Generate an ecological forecast for primary productivity. (Activity A)    
3. Describe how to assess ecological forecast accuracy. (Activity A)     
4. Describe how data assimilation affects forecast accuracy and uncertainty. (Activity B)    
5. Explain how updating models with data collected at different time scales (e.g., daily, weekly) and with different levels of associated uncertainty affects ecological forecasts. (Activity B, C) 
  
## Key Concepts

### What is data assimilation?

Data assimilation is the process of updating models with data. In ecological forecasting, data assimilation is the process of updating ecological forecasting models with new environmental data as they become available.

### How does the amount of uncertainty in model predictions and data affect the process of data assimilation?

The amount of uncertainty in model predictions and data determines how much we adjust our forecasts based on new observations. If we observe a new data point and we have low observation uncertainty (i.e., high confidence in the accuracy of that observation), our forecast starting conditions will be adjusted to closely correspond to the new observation. If we observe a new data point and we have high observation uncertainty (i.e., low confidence in the accuracy of that observation), our forecast starting conditions will not be adjusted as much.

### How does the frequency of observations affect data assimilation?

More frequent observations allow us to update our forecast models more often, potentially improving forecast accuracy.

If you would like to review key slides from the introductory presentation accompanying this module, you can do so by navigating to https://macrosystemseddie.shinyapps.io/module7/ and clicking on the 'Presentation' tab.
  
## Overview

In this module, we will generate one-day-ahead forecasts of lake chlorophyll-a. First, we will generate forecasts that do not assimilate any data. This will involve the following steps:  

#### Activity A
Objective 1. Read in and visualize chlorophyll-a data from Lake Barco, FL, USA.  
Objective 2. Explore autocorrelation of Lake Barco chlorophyll-a data.    
Objective 3. Fit an autoregressive forecast model.   
Objective 4. Generate a one-day-ahead forecast with uncertainty.

Next, we will explore the effect of **data assimilation** on forecast output.  

#### Activity B
Objective 5. Compare one-day-ahead forecasts generated with and without data assimilation.  
Objective 6. Compare one-day-ahead forecasts generated with data assimilation, using data with low vs. high observation uncertainty.  
Objective 7. Compare a series of one-day-ahead forecasts with no data assimilation, weekly data assimilation, and daily data assimilation. 

Then, you will be asked to apply what you have learned about how data collection frequency and observation uncertainty affect data assimilation to improve forecast accuracy.

#### Activity C
Finally, you will have the opportunity to explore the effect of data assimilation on forecasts for a water quality variable of your choice.

##### Independent coding activity
Objective 8. Fit a model and calculate uncertainty for a different water quality variable.  
Objective 9. Determine the optimal frequency of data assimilation for forecasts of your water quality variable.
  
Example code is provided for Objectives 1-7, and you will be asked several short answer questions to interpret code output. Many of these short answer questions parallel (and in some cases are identical to) questions in the R Shiny app version of the module. Questions which are identical to those in the Shiny app will be indicated with **(Shiny)**, while questions unique to this RMarkdown will be indicated with **(Rmd)**. Note that question numbers will differ between the RMarkdown and the Shiny app, even if the question text is the same. Beginning in Objective 8, you will be guided to build on the module code by adjusting example code or developing your own code. Keep in mind that the example code is provided to help you, and use as much of it as you can in completing the questions embedded in Objective 8-9. There are a total of 67 questions. Please see the module rubric for possible points per question and confirm with your instructor whether and how the module will be graded.  

## Feedback

<https://github.com/MacrosystemsEDDIE/module7_R/issues>


## Instructions
  - Attend the introductory PowerPoint lecture provided by your instructor; slides may be downloaded [here](https://d32ogoqmya1dw8.cloudfront.net/files/eddie/teaching_materials/modules/instructors_powerpoint_16626467611382673378.pptx)
  - Open the notebook `assignment/module7.Rmd` in RStudio
  - Work through the exercises described in the notebook.
  - `knit` + commit output files to GitHub

## Context

This module contains code to reproduce the basic functionality of "Macrosystems EDDIE Module 7: Using Data to Improve Ecological Forecasts", found at https://serc.carleton.edu/eddie/teaching_materials/modules/module7.html. The code can be used by students to better understand what is happening "under the hood" of the Module 7 Shiny app, which can be found at the following link:  
https://macrosystemseddie.shinyapps.io/module7/. 
  
Alternatively, students can complete this version of the module instead of the Shiny app version.  

## Timeframe

2 75-minute class periods are allocated to this module

## Background Reading
  
Optional pre-class readings and videos:  
  
**Articles:**  
  
1. Silver, N. 2012. The Signal and the Noise: Why so many Predictions Fail – but some Don't. Penguin Books. 
2. Dietze, M., & Lynch, H. (2019). Forecasting a bright future for ecology. Frontiers in Ecology and the Environment, 17(1), 3. https://doi.org/10.1002/fee.1994. 
3. Dietze, M. C., Fox, A., Beck-Johnson, L. M., Betancourt, J. L., Hooten, M. B., Jarnevich, C. S., Keitt, T. H., Kenney, M. A., Laney, C. M., Larsen, L. G., Loescher, H. W., Lunch, C. K., Pijanowski, B. C., Randerson, J. T., Read, E. K., Tredennick, A. T., Vargas, R., Weathers, K. C., & White, E. P. (2018). Iterative near-term ecological forecasting: Needs, opportunities, and challenges. Proceedings of the National Academy of Sciences, 115(7), 1424–1432. https://doi.org/10.1073/pnas.1710231115 
  
**Videos:**

1. NEON's [Ecological Forecast: The Science of Predicting Ecosystems](https://www.youtube.com/watch?v=Lgi_e7N-C8E&t=196s&pbjreload=101)  
2. Fundamentals of Ecological Forecasting Series
      - [Why Forecast?](https://www.youtube.com/watch?v=kq0DTcotpA0&list=PLLWiknuNGd50Lc3rft4kFPc_oxAhiQ-6s&index=1)
      - [Uncertainty Analysis](https://www.youtube.com/watch?v=rDCkjzVQNSw&list=PLLWiknuNGd50Lc3rft4kFPc_oxAhiQ-6s&index=12)
  
## References
  
This module is derived from Lofton, M.E., Moore, T. N., Carey, C.C. and Thomas, R. Q. 07 March 2024. Macrosystems EDDIE: Using Data to Improve Ecological Forecasts. Macrosystems EDDIE Module 7, Version 1. http://module7.macrosystemseddie.org. Module development was supported by NSF grants DEB-1926050 and DBI-1933016.  
  
-   author: Mary Lofton (@melofton)
-   contact: [melofton\@vt.edu](mailto:melofton@vt.edu)
-   url: https://serc.carleton.edu/eddie/teaching_materials/modules/module6.html
-   date: 2023-12-12
-   license: MIT, CC-BY
-   copyright: Mary Lofton

