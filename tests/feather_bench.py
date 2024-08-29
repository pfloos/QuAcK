
from molecule import Molecule


He = Molecule(
    name="He",
    multiplicity=1,
    geometry=[
        {"element": "He", "x": 0.0, "y": 0.0, "z": 0.0}
    ],
    properties={
        "properties_rhf":{
            "6-31g": {
                "RHF energy": -2.855160426884076,
                "RHF HOMO energy": -0.914126628614305,
                "RHF LUMO energy": 1.399859335225087,
                "RHF dipole moment": 0.000000000000000,
                "RMP2 correlation energy": -0.011200122910187,
                "CCD correlation energy": -0.014985063408247,
                "DCD correlation energy": -0.014985062907429,
                "CCSD correlation energy": -0.015001711549550,
                "drCCD correlation energy": -0.018845374502248,
                "rCCD correlation energy": -0.016836324636164,
                "crCCD correlation energy": 0.008524677369855,
                "lCCD correlation energy": -0.008082420815100,
                "pCCD correlation energy": -0.014985062519068,
                "RCIS singlet excitation energy": 1.911193619935257,
                "RCIS triplet excitation energy": 1.455852629402236,
                "phRRPA correlation energy": -0.018845374129105,
                "phRRPAx correlation energy": -0.015760565121283,
                "crRRPA correlation energy": -0.008868581132405,
                "ppRRPA correlation energy": -0.008082420815100,
                "RG0F2 correlation energy": -0.011438430540374,
                "RG0F2 HOMO energy": -0.882696116247871,
                "RG0F2 LUMO energy": 1.383080391811630,
                "evRGF2 correlation energy": -0.011448483158486,
                "evRGF2 HOMO energy": -0.881327878713477,
                "evRGF2 LUMO energy": 1.382458968133448,
                "RG0W0 correlation energy": -0.019314094399756,
                "RG0W0 HOMO energy": -0.870533880190454,
                "RG0W0 LUMO energy": 1.377171287010956,
                "evRGW correlation energy": -0.019335511771724,
                "evRGW HOMO energy": -0.868460640957913,
                "evRGW LUMO energy": 1.376287581471769,
                "RG0T0pp correlation energy": -0.008082420815100,
                "RG0T0pp HOMO energy": -0.914126628614305,
                "RG0T0pp LUMO energy": 1.399859335225087,
                "evRGTpp correlation energy": -0.008082420815100,
                "evRGTpp HOMO energy": -0.914126628614305,
                "evRGTpp LUMO energy": 1.399859335225087
            }
        },
        "properties_uhf":{
            "6-31g": {
            }
        },
        "properties_ghf":{
            "6-31g": {
            }
        },
        "properties_rohf":{
            "6-31g": {
            }
        }
    }
)

# ---

#H2O = Molecule(
#    name="H2O",
#    multiplicity=1,
#    geometry=[
#        {"element": "O", "x":  0.0000, "y": 0.0000, "z": 0.0000},
#        {"element": "H", "x":  0.7571, "y": 0.0000, "z": 0.5861},
#        {"element": "H", "x": -0.7571, "y": 0.0000, "z": 0.5861}
#    ],
#    properties={
#        "cc-pvdz": {
#    }
#)

# ---

FeatherBench = [
    He, 
    #H2O
]





