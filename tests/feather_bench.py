
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
                "RHF energy": -2.855160426154444,
                "RHF HOMO energy": -0.914126628640145,
                "RHF LUMO energy": 1.399859335255765,
                "RHF dipole moment": 0.0,
                "MP2 correlation energy": -0.011200122909934,
                "CCD correlation energy": -0.014985063116,
                "CCSD correlation energy": -0.015001711549092,
                "drCCD correlation energy": -0.01884537385338,
                "rCCD correlation energy": -0.016836322809386,
                "crCCD correlation energy": 0.008524676641474,
                "lCCD correlation energy": -0.00808242082105,
                "CIS singlet excitation energy": 1.911193619991987,
                "CIS triplet excitation energy": 1.455852629458543,
                "phRPA correlation energy": -0.018845374128748,
                "phRPAx correlation energy": -0.015760565120758,
                "crRPA correlation energy": -0.008868581132249,
                "ppRPA correlation energy": -0.008082420814972,
                "G0F2 correlation energy": -0.011438430540104,
                "G0F2 HOMO energy": -0.882696116274599,
                "G0F2 LUMO energy": 1.383080391842522,
                "G0W0 correlation energy": -0.019314094399372,
                "G0W0 HOMO energy": -0.87053388021722,
                "G0W0 LUMO energy": 1.377171287041735,
                "evGW correlation energy": -0.019335511771337,
                "evGW HOMO energy": -0.868460640984803,
                "evGW LUMO energy": 1.376287581502582,
                "G0T0pp correlation energy": -0.008161908540634,
                "G0T0pp HOMO energy": -0.898869172597701,
                "G0T0pp LUMO energy": 1.383928087417952,
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
    }
)

# ---

H2O = Molecule(
    name="H2O",
    multiplicity=1,
    geometry=[
        {"element": "O", "x":  0.0000, "y": 0.0000, "z": 0.0000},
        {"element": "H", "x":  0.7571, "y": 0.0000, "z": 0.5861},
        {"element": "H", "x": -0.7571, "y": 0.0000, "z": 0.5861}
    ],
    properties={
        "properties_rhf":{
            "cc-pvdz": {
                "RHF energy": -85.21935817501823,
                "RHF HOMO energy": -0.493132793449897,
                "RHF LUMO energy": 0.185534869842355,
                "RHF dipole moment": 0.233813698748474,
                "MP2 correlation energy": -0.203978216774657,
                "CCD correlation energy": -0.212571260121257,
                "CCSD correlation energy": -0.213302190845899,
                "drCCD correlation energy": -0.231281853419338,
                "rCCD correlation energy": -0.277238348710547,
                "crCCD correlation energy": 0.18014617422324,
                "lCCD correlation energy": -0.15128653432796,
                "CIS singlet excitation energy": 0.338828950934568,
                "CIS triplet excitation energy": 0.304873339484139,
                "phRPA correlation energy": -0.231281866582435,
                "phRPAx correlation energy": -0.310796738307943,
                "crRPA correlation energy": -0.246289801609294,
                "ppRPA correlation energy": -0.151286536255888,
                "G0F2 correlation energy": -0.217807591229668,
                "G0F2 HOMO energy": -0.404541451101377,
                "G0F2 LUMO energy": 0.16650398400197,
                "G0W0 correlation energy": -0.23853664665404,
                "G0W0 HOMO energy": -0.446828623007469,
                "G0W0 LUMO energy": 0.173026609033024,
                "evGW correlation energy": -0.239414217281308,
                "evGW HOMO energy": -0.443076613314424,
                "evGW LUMO energy": 0.172691758111392,
                "G0T0pp correlation energy": -0.156214864467344,
                "G0T0pp HOMO energy": -0.452117482732615,
                "G0T0pp LUMO energy": 0.16679206983464,
            }
        }
    }
)

# ---

FeatherBench = [
    He, 
    H2O
]





