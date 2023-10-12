#TNDDR: Efficient and doubly robust estimation of COVID-19 vaccine effectiveness under the test-negative design

This repository is the official implementation of [Efficient and doubly robust estimation of COVID-19 vaccine effectiveness under the test-negative design](https://arxiv.org/abs/2310.04578).

## Introduction

This study focuses on the derivation and evaluation of an efficient nonparametric estimator of vaccine effectiveness under an outcome-dependent sampling design known as the testnegative design (TND). This relatively new epidemiological design is currently being used globally to evaluate vaccine effectiveness for COVID, seasonal influenza, and other infectious diseases. It is extremely useful in the context of rapid vaccine evaluation since
it can be deployed in essentially real time. 

We derived the efficiency theory under the distribution arising from the biased TND sampling and a resulting one-step, doubly robust, and locally efficient estimator that can be implemented with machine learning methods. The properties of the proposed estimator are demonstrated theoretically and through simulation study. We then employ it to estimate the effectiveness of COVID-19 booster vaccination in a dataset of healthcare workers from the province of Québec, Canada.

## Requirements

 - R 3.6
 - `earth`
 - `sandwich`
 - `hal9001`
 - `dplyr`
 - `haven`
 - `sas7bdat`
 
## Contents

  1. `README.txt`: implementation details of source code and contents 

  2. Source codes of TNDDR and data generation environment
  
     a). `TureValue.R`: Calculate the true mRR values for different settings.
     
     Folder Study 1: codes for Study 1 (assessing the efficiency of TNDDR)

     1). `Study1.R`: funcitons for Table 1 that presents the results for mRR estimation and inference; 

     2). `Study1RootN.R`: main functions for Figure 2, which illustrates the average error of estimators across sample sizes;
     
     Folder Study 2: codes for Study 2 (validating the double robustness of TNDDR)
     
     1). `Study2Table2.R`: codes for Table 2, which displays results from the four scenarios, with different sample sizes 
     
     Folder Realdata: codes for Realdata analysis
     
 

