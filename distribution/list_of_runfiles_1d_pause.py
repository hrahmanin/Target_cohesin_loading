import glob
import ast 
import datetime

LEF_pause= 0.9

pause_multiplier = 1/(1-LEF_pause)

lifetimes_b = [66]
lifetimes = [pause_multiplier*life for life in lifetimes_b]

velocities = [1]

face_stalls = [1.0]
back_stalls = [0.0]

CTCF_lifetimes_b =  [ 30]
CTCF_lifetimes = [pause_multiplier*clife for clife in CTCF_lifetimes_b]
CTCF_offtimes_b = [18]
CTCF_offtimes = [pause_multiplier*cof for cof in CTCF_offtimes_b]

stall_dists = [100]
LEF_birth = 0.1

LEF_separations = [75]
sites_per_monomer = 10
replication_number = 10
monomer_per_replica = 1000
occup_list = [0.1]
steps_b = 200
steps = steps_b/pause_multiplier

already_processed = []
for fname  in glob.glob('/home1/start/polychrom/projects/Dynamic_boundary_elements/simulations/sims/*_STALL*'):
    already_processed.append(fname.split('/')[-1])

with open(str(datetime.date.today())+'_runfile.txt','w') as f:
    for lifetime in lifetimes:
        for velocity_multiplier in velocities:
            for face_stall in face_stalls:
                for back_stall in back_stalls:
                    for LEF_separation in LEF_separations:
                        for CTCF_lifetime in CTCF_lifetimes:
                            for CTCF_offtime in CTCF_offtimes:
                                paramset = (
                                    f'folder_face_{face_stall:.1f}'
                                    f'_back_{back_stall:.1f}'
                                    f'_Clife_{CTCF_lifetime:.0f}'
                                    f'_Cof_{CTCF_offtime:.0f}'
                                    f'_life_{lifetime:.1f}'
                                    f'_slife_{lifetime:.1f}'
                                    f'_birth_{LEF_birth:.1f}'
                                    f'_pause_{LEF_pause:.1f}'
                                    f'_sep_{LEF_separation:.0f}'
                                    f'_site_{sites_per_monomer:.0f}'
                                    f'_replica_{replication_number:.0f}'
                                    f'_vel_{velocity_multiplier:.0f}'
                                    f'_steps_{steps:.0f}'
                                )
                                if paramset not in already_processed:
                                    f.write(paramset + '\n')
                                else:
                                    print('already done')

                                
            

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }
