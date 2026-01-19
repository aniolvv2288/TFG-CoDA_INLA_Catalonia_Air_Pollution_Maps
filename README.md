# Final Degree Project (TFG) – Compositional Data Analysis (CoDA) for Air Pollution Mapping in Catalonia

This repository contains the code, data workflows, and documentation for my Final Degree Project (TFG) in Applied Statistics at Universitat Autònoma de Barcelona (UAB).

## Project Overview

The main goal of this project is to characterize and map atmospheric pollution across Catalunya using Compositional Data Analysis (CoDA). The work produces continuous environmental exposure surfaces that not only describe the spatial distribution of pollutants but also provide a robust foundation for future epidemiological research, including studies on the relationship between air pollution exposure and mental health or cardiovascular outcomes.

## Key Components

- **Data Cleaning & Preprocessing**  
  Scripts for importing, cleaning, formatting, and preparing environmental datasets.

- **Data Imputation**  
  Missing values handled using KNN imputatio from the *robCompositions* package, selected after a sensitivity comparison with robust LTS regression.

- **Outliers Exploration**  
  Outlier detection performed with *OutlierClassifier1* from the *compositions* package, following the compositional outlier‑detection framework of Tolosana‑Delgado & Mueller (2021).

- **Compositional Data Analysis**  
  Application of CoDA transformations (clr, ilr), balances, and multivariate modelling to analyse pollutant compositions.

- **Exploratory Data Analysis**  
  Exploratory assessment of pollutant compositions using the variation matrix, multivariate PCA on clr‑transformed data, and T‑space interpretation.

- **Exposure Mapping**  
  Construction of continuous exposure surfaces using a Bayesian hierarchical spatial model implemented with INLA + SPDE to approximate Gaussian random fields.

- **Reproducibility**  
  Modular R scripts, version control via GitHub, and documentation ensuring a transparent and reproducible workflow.

## Structure

- `/data`: Raw and processed datasets (with `.gitignore` for sensitive files)
- `/scripts`: R scripts for each analysis step
- `/output`: Maps, tables, and results

## Authors

- Author: Aniol Vilamala Vidal
- Tutor: David Moriña Soler


Universitat Autònoma de Barcelona (UAB)
