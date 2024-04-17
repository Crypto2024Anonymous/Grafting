# Grafting
This is an anonymous git repository for the Crypto 2024 submission #192 "Grafting: Complementing RNS in CKKS."

The codes are originally from [Lattigo](https://github.com/tuneinsight/lattigo) v5.0.2. commit `4cce9a48c1daaa2dd122921822f5ad70cd444156`.

Try 
```
go run compare.go
```
for comparing the run-time of **Original Lattigo Bootstrapping** versus **Out-sourced Bootstrapping based on Grafting**. 

The **Out-sourced Bootstrapping based on Grafting** is an extension of the linear transformation in the word-sized moduli chain in the submission paper. 
Currently, it does not include the last modulus switching (which reported 0.3% of the total bootstrapping time in HEaaN implementation). 
A complete version will be uploaded in few days. 
