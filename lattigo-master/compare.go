// Package main implements an example showcasing the comparison of the two CKKS bootstrappings:
// 1) original bootstrapping implemented in Lattigo with example 128-bit secure parameters logN=16, max logPQ=638, and 1520 for ciphertexts and keys, respectively,
// 2) bootstrapping with Grafting (a.k.a. Out-sourced BTS) in the word-size moduli chain, but without the last modulus switching back to the original moduli chain (under developing)
// Note: the last modulus switching recorded a run-time of 0.3% of the total BTS, which is negligible for this comparison.
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16
	var flagShortInt int
	flagShortInt = 0

	if *flagShort {
		flagShortInt = 1
		LogN -= 3
	}

	//===========================================
	//=== 1) RESIDUAL PARAMETERS for two BTSs ===
	//===========================================

	// Define the residual parameters, used outside of the bootstrapping circuit.
	// LogN=16, logQ = 55 + 10*40 = 55 + 6*60 + 40 and logP = 3*61, so LogQP = 638, well over 128-bit of security.

	// Original parameters:
	OriginParams, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,                                              // Log2 of the ring degree
		LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61},                                 // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 40,                                                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})
	if err != nil {
		panic(err)
	}

	// Out-sourced BTS parameters:
	OutsourcedParams, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,                                  // Log2 of the ring degree
		LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 40}, // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61},                     // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 40,                                    // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})
	if err != nil {
		panic(err)
	}

	fmt.Println("*****************************************************************")
	fmt.Println("* Comparison of the two RNS-CKKS bootstrappings:                *")
	fmt.Println("* 1) Original Lattigo BTS,                                      *")
	fmt.Println("* 2) Out-sourced BTS following Grafting paper                   *")
	fmt.Println("* Note, the last ModSwitch (0.3%% of BTS runi-time) is omitted. *")
	fmt.Println("*****************************************************************\n")

	fmt.Println("Original BTS:")
	timeElapedOriginal := BenchBTS(OriginParams, LogN, flagShortInt)
	fmt.Println("Out-sourced BTS:")
	timeElapedOutsourced := BenchBTS(OutsourcedParams, LogN, flagShortInt)

	fmt.Printf("*********************\n* Comparison Result *\n*********************\n\n")
	fmt.Printf("Original vs. Grafting BTS timings: %.2f ms vs. %.2f ms", timeElapedOriginal, timeElapedOutsourced)
	fmt.Printf("\nOut-sourced BTS is %fx faster!!!\n", timeElapedOriginal/timeElapedOutsourced)
	fmt.Printf("Note, the Out-sourced BTS do not include the last ModSwitch, which is around 0.3%% of the total run-time.\n\n")
}

func BenchBTS(params hefloat.Parameters, LogN int, flagShortInt int) (timeElapsed float64) {

	//=== 2) BOOTSTRAPPING PARAMETERSLITERAL ===
	btpParametersLit := bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(LogN),
		LogP: []int{61, 61, 61, 61},
		Xs:   params.Xs(),
	}

	//=== 3) BOOTSTRAPPING PARAMETERS ===
	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if flagShortInt == 0 {
		// Corrects the message ratio Q0/|m(X)| to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	// We print some information about the residual parameters.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(), params.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	// And some information about the bootstrapping parameters.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//=== 4) KEYGEN & ENCRYPT ===
	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	encoder := hefloat.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//=== 5) BOOTSTRAPPING ===
	// Instantiates the bootstrapper
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	// Generate a random plaintext with values uniformely distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, params.MaxSlots())
	for i := range valuesWant {
		valuesWant[i] = sampling.RandComplex128(-1, 1)
	}

	// We encrypt at level 0
	plaintext := hefloat.NewPlaintext(params, 0)
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext1, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	// fmt.Println()
	// fmt.Println("*Precision of values vs. ciphertext")
	// valuesTest1 := printDebug(params, ciphertext1, valuesWant, decryptor, encoder)
	valuesTest1 := make([]complex128, ciphertext1.Slots())
	if err := encoder.Decode(decryptor.DecryptNew(ciphertext1), valuesTest1); err != nil {
		panic(err)
	}

	// Bootstrap
	fmt.Println("Bootstrapping...")
	start := time.Now()
	ciphertext2, err := eval.Bootstrap(ciphertext1)
	if err != nil {
		panic(err)
	}
	elapsed := time.Since(start)
	fmt.Printf("Done! BTS took %s\n\n", elapsed)

	timeElapsed = float64(elapsed.Milliseconds())

	//==================
	//=== 6) DECRYPT ===
	//==================

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("*Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext2, valuesTest1, decryptor, encoder)

	return
}

func printDebug(params hefloat.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *hefloat.Encoder) (valuesTest []complex128) {

	valuesTest = make([]complex128, ciphertext.Slots())

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := hefloat.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())
	// fmt.Println()

	return
}
