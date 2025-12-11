from ase.calculators.orca import ORCA
from ase.calculators.orca import OrcaProfile


def ORCA_calc():
    MyOrcaProfile = OrcaProfile(command='/software/ORCA/6.0.0/icelake/orca')

    calc = ORCA(profile=MyOrcaProfile,
        charge=0,
        mult=1, 
        orcasimpleinput='UKS B3LYP 6-31G* D3ZERO ENGRAD PAL8',
        orcablocks="""
    %scf
    Guess MORead
    MOInp "/Path/to/brokensym/orbitals/guess.gbw"
    end
    """ 
    )
    return calc
