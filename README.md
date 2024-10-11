## Testing the Permanent Income Hypothesis (PIH) Using Consumption and Income Data

### Intro
As part of my undergraduate degree, I studied the Quantitative Economics (QE) module which is essentially an introductory econometrics course that follows the renowned textbook by Stock and Watson. This meant that I had the opportunity to build on my pre-existing foundations in R that I developed during Q-Step 1 and 2 from the Politics Department. While we did do many excercises within R for QE, this was limited almost solely to R fundamentals and then regression analysis using Ordinary Least Squares and Instrumental Variables, which left the time series portion of the module without any practical experience in its application. While we did go into depth about the theory and interpretation of results and tutors suggested packages that we could use to carry out a hypothetical project, I felt that it would be beneficial to build my own practical understanding. So, in my free time, I took a look at different online resources and built a project in RStudio to do just that. At the time, I was revising Macroeconomics while preparing for the Money & Banking module, which is why I thought that testing the implications of the Permanent Income Hypothesis would be perfect as the subject of my analysis.  
Although the investigation is not very rigorous nor ground-breaking, I feel that I achieved my goal of learning the basics of time series analysis in R, while also gaining some firsthand insight into how well the Permanent Income Hypothesis aligns with Macroeconomic data.

### Abstract
The Permanent Income Hypothesis (PIH), proposed by Milton Friedman, suggests that individuals tend to smooth their consumption over time based on their expected "permanent" income, rather than on temporary fluctuations in currrent income. In this way, temporary changes in income should not lead to significant changes in consumption, so long as their long-term income expectations (on which their spending decisions are based) remain stable. This project tests the PIH by exploring the relationship between consumption and income using time series analysis techniques on aggregate-level data from FRED.

### Outline

Our main objectives will be:
- Investigate the dynamics between income and consumption, testing for Granger causality
- Explore structural breaks, especially around economic crises or major policy changes
- Test for a cointegrating relationship between income and consumption
- Forecast future consumption and income

Key steps and techniques:
1. Data collection

2. Exploratory data analysis
    - Visualise the time series
    - Plot ACFs to detect persistence and seasonality

3. Transform to stationarity
    - Use log-transformed data to stabilise the variance and consider growth-rates
    - Determine order of integration using ADF tests

4. Model selection
    - Select lag lengths using information criteria for AR and ADL models
    
5. Test for Granger causality

6. Test for structural breaks using QLR test

7. Re-evaluate models and forecast

8. Test for cointegration
    - Estimate cointegrating relationship
    - Test for stationarity in residuals
