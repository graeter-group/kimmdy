#
#  1  Bond             2  Morse            3  Angle            4  Proper-Dih.   
#  5  Per.-Imp.-Dih.   6  LJ-14            7  Coulomb-14       8  LJ-(SR)       
#  9  Disper.-corr.   10  Coulomb-(SR)    11  Coul.-recip.    12  Potential     
# 13  Kinetic-En.     14  Total-Energy    15  Conserved-En.   16  Temperature   
# 17  Pres.-DC        18  Pressure        19  dVremain/dl     20  Box-X         
# 21  Box-Y           22  Box-Z           23  Volume          24  Density       
# 25  pV              26  Enthalpy        27  Vir-XX          28  Vir-XY        
# 29  Vir-XZ          30  Vir-YX          31  Vir-YY          32  Vir-YZ        
# 33  Vir-ZX          34  Vir-ZY          35  Vir-ZZ          36  Pres-XX       
# 37  Pres-XY         38  Pres-XZ         39  Pres-YX         40  Pres-YY       
# 41  Pres-YZ         42  Pres-ZX         43  Pres-ZY         44  Pres-ZZ       
# 45  #Surf*SurfTen   46  Box-Vel-XX      47  Box-Vel-YY      48  Box-Vel-ZZ    
# 49  T-Protein                           50  T-non-Protein                     
# 51  Lamb-Protein                        52  Lamb-non-Protein                  

# %%
# with pairs, PR current state
# %%
kimmdy

# %%
gmx energy -f ./alanine_hat_000/5_relax/relax.edr

# %%
kimmdy-analysis energy alanine_hat_v3 -o --terms potential lj-14 'LJ-(SR)' coulomb-14 'coulomb-(sr)' 'coul.-recip'

# %%
kimmdy --input ./kimmdy-simplified.yml

# %%
kimmdy-analysis energy alanine_hat_v1_emulated_by_v3 -o --terms potential lj-14 'LJ-(SR)' coulomb-14 'coulomb-(sr)' 'coul.-recip'

# %%
diff ./alanine_hat_001/9_apply_recipe/Ala_out_before.top ./alanine_hat_001/9_apply_recipe/Ala_out_after.top
# ./alanine_hat_001/9_apply_recipe/Ala_out_relax.top
diff ./alanine_hat_001/9_apply_recipe/Ala_out_before_with_solvent_bonds.top ./alanine_hat_001/9_apply_recipe/Ala_out_before.top

# %%
kimmdy-analysis energy alanine_hat_v1 -o --terms potential lj-14 'LJ-(SR)' coulomb-14 'coulomb-(sr)' 'coul.-recip'

# %%
kimmdy-analysis energy alanine_hat_v2 -o --terms potential lj-14 'LJ-(SR)' coulomb-14 'coulomb-(sr)' 'coul.-recip'

# %%
rm -r ./alanine_hat_00*
