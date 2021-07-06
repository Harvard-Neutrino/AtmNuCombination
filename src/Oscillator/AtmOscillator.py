import nuSQuIDS as nsq
import nuSQUIDSTools


def SMOsc(nuSQ,flavor,in_state,th13,th23,th12,dm21,dm32,dcp,mh):
# Mixing angles in radians
	nuSQ.Set_MixingAngle(0,1,th12)
	nuSQ.Set_MixingAngle(0,2,th13)
	nuSQ.Set_MixingAngle(1,2,th23)
	nuSQ.Set_SquareMassDifference(1,dm21)

	if mh==1:
		nuSQ.Set_SquareMassDifference(2,dm32)
	elif mh==-1:
		nuSQ.Set_SquareMassDifference(2,mh*dm32+dm21)
	else:
		print('No valid MH value, setting MH to Normal')
		nuSQ.Set_SquareMassDifference(2,dm32)

	nuSQ.Set_CPPhase(0,2,dcp)

	nuSQ.Set_initial_state(in_state,nsq.Basis.flavor)
	nuSQ.EvolveState()

	j = int(abs(flavor) / 2) % 6

	return nuSQ.EvalFlavor(j)
