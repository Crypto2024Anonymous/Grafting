package ckks

import (
	// "math/bits"
	"github.com/tuneinsight/lattigo/v5/ring"
	// "github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	// "github.com/tuneinsight/lattigo/v5/utils"
)

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type ModulusSwitcher struct {
	*ModulusSwitcherBase
	*ModulusSwitcherBuffers

	PermuteNTTIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type ModulusSwitcherBase struct {
	params1 Parameters
	params2 Parameters
}

type ModulusSwitcherBuffers struct {
	BuffQ1 [2]ring.Poly
	BuffQ2 [2]ring.Poly
}

func newModulusSwitcherBase(params1, params2 Parameters) *ModulusSwitcherBase {
	rs := new(ModulusSwitcherBase)
	rs.params1 = params1
	rs.params2 = params2
	return rs
}

func newModulusSwitcherBuffers(params1, params2 Parameters) *ModulusSwitcherBuffers {

	buff := new(ModulusSwitcherBuffers)
	buff.BuffQ1 = [2]ring.Poly{params1.RingQ().NewPoly(), params1.RingQ().NewPoly()}
	buff.BuffQ2 = [2]ring.Poly{params2.RingQ().NewPoly(), params2.RingQ().NewPoly()}

	return buff
}

func NewModulusSwitcher(params1, params2 Parameters) (rs *ModulusSwitcher) {
	rs = new(ModulusSwitcher)

	rs.ModulusSwitcherBase = newModulusSwitcherBase(params1, params2)
	rs.ModulusSwitcherBuffers = newModulusSwitcherBuffers(params1, params2)

	// Warning! omit checking RingP of srcparams and dstparams are the same

	if params1.RingQ() != nil && params2.RingQ() != nil {
		rs.BasisExtender = ring.NewBasisExtender(params1.RingQ(), params2.RingQ())
		rs.Decomposer = ring.NewDecomposer(params1.RingQ(), params2.RingQ())
	}

	return
}
func (rs *ModulusSwitcher) SwitchSecretKey(kgenIn, kgenOut *rlwe.KeyGenerator, skIn *rlwe.SecretKey) (skOut *rlwe.SecretKey) {
	var ringQ1 *ring.Ring
	var ringQ2 *ring.Ring
	var buff *[2]ring.Poly
	// var params_src Parameters
	var params_dst Parameters

	if skIn.LevelQ() == rs.params1.MaxLevelQ() {
		// params_src = rs.params1
		params_dst = rs.params2
		ringQ1 = rs.params1.RingQ()
		ringQ2 = rs.params2.RingQ()
		buff = &rs.BuffQ1

	} else if skIn.LevelQ() == rs.params2.MaxLevelQ() {
		// params_src = rs.params2
		params_dst = rs.params1
		ringQ1 = rs.params2.RingQ()
		ringQ2 = rs.params1.RingQ()
		buff = &rs.BuffQ2

	} else {
		panic("skIn LevelQ does not match any of the MaxLevelQ of the parameters")
	}

	skOut = rlwe.NewSecretKey(params_dst)
	rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ1, ringQ2, skIn.Value.Q, buff[0], skOut.Value.Q)
	rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ1, rs.params1.RingP(), skIn.Value.Q, buff[0], skOut.Value.P)
	return skOut
}

func (rs *ModulusSwitcher) SwitchNew(ctIn *rlwe.Ciphertext, dir bool) (ciphertextOut *rlwe.Ciphertext) {

	var ringQ1 *ring.Ring
	var ringQ2 *ring.Ring
	var levelQ2 int
	var buffQ1 *[2]ring.Poly
	var buffQ2 *[2]ring.Poly
	var ctOut *rlwe.Ciphertext
	levelQ1 := ctIn.Level()

	if dir {
		buffQ1 = &rs.BuffQ1
		buffQ2 = &rs.BuffQ2
		levelQ2 = levelQ1 - rs.params1.RingQ().Level() + rs.params2.RingQ().Level()
		ctOut = rlwe.NewCiphertext(rs.params1, ctIn.Degree(), levelQ1)
		ringQ1 = rs.params1.RingQ().AtLevel(levelQ1)
		ringQ2 = rs.params2.RingQ().AtLevel(levelQ2)
		ctOut.Scale = rs.params2.DefaultScale()
		// } else if ctIn.LogScale() == float64(rs.params2.LogDefaultScale()) {
	} else {
		buffQ1 = &rs.BuffQ2
		buffQ2 = &rs.BuffQ1
		levelQ2 = levelQ1 - rs.params2.RingQ().Level() + rs.params1.RingQ().Level()
		ctOut = rlwe.NewCiphertext(rs.params1, ctIn.Degree(), levelQ2)
		ringQ1 = rs.params2.RingQ().AtLevel(levelQ1)
		ringQ2 = rs.params1.RingQ().AtLevel(levelQ2)
		ctOut.Scale = rs.params1.DefaultScale()
	}

	ctOut.El().IsNTT = ctIn.El().IsNTT
	ctOut.El().IsBatched = ctIn.El().IsBatched
	ctOut.El().LogDimensions = ctIn.El().LogDimensions

	if dir {
		for i := range ctIn.Value {
			ringQ1.INTTLazy(ctIn.Value[i], buffQ1[i])
			rs.BasisExtender.ModUpQtoP(levelQ1, levelQ2, buffQ1[i], buffQ2[i])
			ringQ2.NTT(buffQ2[i], ctOut.Value[i])
		}
	} else {
		for i := range ctIn.Value {
			ringQ1.INTTLazy(ctIn.Value[i], buffQ1[i])
			rs.BasisExtender.ModUpPtoQ(levelQ1, levelQ2, buffQ1[i], buffQ2[i])
			ringQ2.NTTLazy(buffQ2[i], ctOut.Value[i])
		}
	}
	return ctOut
}
