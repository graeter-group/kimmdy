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
kimmdy --input kimmdy.yml

# %%
kimmdy-analysis energy alanine_hat_000 -o --terms lj-14 'lj-(sr)' total-energy

# %%
kimmdy-analysis energy alanine_hat_000 -o --terms bond angle proper-dih

# %%
kimmdy-analysis energy alanine_hat_000 -o --terms potential kinetic-en total-energy

# %%
kimmdy-analysis energy alanine_hat_016 -o --terms temperature potential

# %%
kimmdy-analysis energy alanine_hat_018 -o --terms temperature potential

# %%
kimmdy-analysis trjcat alanine_hat_000 -o -g0

# %%
echo 'lj-sr:protein-protein' | gmx energy -f ./alanine_hat_000/1_equilibrium/equilibrium.edr -o energy1.xvg
echo 'lj-sr:protein-protein' | gmx energy -f ./alanine_hat_000/5_relax/relax.edr -o energy2.xvg
echo 'lj-sr:protein-protein' | gmx energy -f ./alanine_hat_000/6_equilibrium/equilibrium.edr -o energy3.xvg
xmgrace ./energy*.xvg

#%%
echo 'lj-14:protein-protein' | gmx energy -f ./alanine_hat_000/1_equilibrium/equilibrium.edr -o energy1.xvg
echo 'lj-14:protein-protein' | gmx energy -f ./alanine_hat_000/5_relax/relax.edr -o energy2.xvg
echo 'lj-14:protein-protein' | gmx energy -f ./alanine_hat_000/6_equilibrium/equilibrium.edr -o energy3.xvg

xmgrace ./energy*.xvg
