# BDT Tree Model for Mortgage-Backed Securities

This project applies the **Black-Derman-Toy (BDT)** short-rate interest rate model and **Nelson-Siegel** term structure to evaluate the pricing and risk characteristics of mortgage-backed securities (MBS), accounting for borrower refinancing behavior.

---

## Project Overview

- **Language**: MATLAB  
-  **Files**:
  - `BDT_MBS.m`: Main script that implements BDT calibration, yield fitting, and MBS valuation
  - `yieldAndVolatility2000.csv`: Input term structure data
-  **Models Used**:
  - Black-Derman-Toy (BDT) short-rate interest rate model
  - Nelson-Siegel yield curve estimation

---

## 🔬 Part I – Building the Interest Rate Framework

###  1. Replicating the BDT Tree
- Reproduced a 3-step BDT tree using `fminsearch` to calibrate short rates and volatility
- Constructed interest rate tree by exponential scaling of volatility across nodes
- Discounted bond cash flows using the `BondTree` function

###  2. Nelson-Siegel Yield and Volatility Curve Fitting
- Estimated yield curve using the Nelson-Siegel functional form
- Calibrated parameters (θ₀, θ₁, θ₂, λ) by minimizing squared error to observed yields
- Extended the term structure up to 360 months

---

##  Part II – MBS Valuation and Risk Analysis

###  1. Monthly Payment and Principal
- Calculated mortgage monthly payments using standard annuity formulas
- Tracked outstanding principal over time

###  2. MBS Valuation: With vs Without Prepayment

####  With Prepayment:
- Modeled borrower behavior to refinance when interest rates fall
- At each node, checked if refinancing was optimal (principal < continuation value)
- Present value of prepayable MBS: **$90,631,628.70**

#### Without Prepayment:
- Assumed fixed cash flows over 30 years
- Valued MBS by backward discounting using the BDT tree
- Present value: **$102,812,530.30**

---

##  Duration and Convexity Analysis

###  Duration
- Measures price sensitivity to interest rate changes
- Without prepayment: duration decreases smoothly with rising rates
- With prepayment: duration drops sharply at low rates due to early refinancing

### 🔹 Convexity
- Captures second-order sensitivity to rates
- Non-prepayable MBS: positive convexity throughout
- Prepayable MBS: **negative convexity** at low rates, flattening upside potential

---

## Discussion on Refinancing Risk

- **Pricing Impact**:
  - Prepayment lowers expected cash flow stability
  - In low-rate environments, prepayable MBS underperform due to refinancing
- **Bank Risk Management**:
  - Hedging via swaps/options
  - Prepayment penalties
  - Upfront fees to price in borrower optionality

---

##  Key Takeaways

- The BDT model, when combined with Nelson-Siegel curves, effectively captures yield dynamics
- Prepayment introduces significant optionality, which lowers MBS valuation and alters risk profiles
- Duration and convexity behave asymmetrically under refinancing assumptions

---

##  Authors

- **Matthieu Fouilloux**  
- **Giulia Gambaretto**  
McGill University 

