# Grafting
This is an anonymous git repository for the Crypto 2024 submission #192 "Grafting: Complementing RNS in CKKS."

The codes are originally from [Lattigo](https://github.com/tuneinsight/lattigo) v5.0.2. commit `4cce9a48c1daaa2dd122921822f5ad70cd444156`.

Try 
```
go run Outsourced_BTS_ModRaise_first.go
```
and
```
go run Outsourced_BTS_Slot2Coeff_first.go
```
for comparing the run-time of **Original Lattigo Bootstrapping** versus **Outsourced Bootstrapping based on Grafting** in two different CKKS bootstrapping settings:
1) ModRaise - Coeff2Slot - EvalMod - Slot2Coeff,
2) Slot2Coeff - ModRaise - Coeff2Slot - EvalMod.

Each reported 1.33x and x1.32 better bootstrapping performance, respectively, averaged from 10 executions (1.4 GHz, Intel Core i5, single core). 

The **Outsourced Bootstrapping based on Grafting** is an extension of the linear transformation in the word-sized moduli chain in the submission paper. 
Currently, 'Outsourced_BTS_ModRaise_first.go' does not include the last modulus switching (which reported 0.3% of the total bootstrapping time in HEaaN implementation), however, 
'Outsourced_BTS_Slot2Coeff_first.go' includes it and provides a squaring and decrypt. 
