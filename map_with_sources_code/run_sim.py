def runSimulation(sim, obs, events:list, seed:int) -> tuple:
    '''
    Function to run CRPropa simulation for H, He, N, Fe with predefined params:
    - 10% energy uncertainty
    - 1 degree directional uncertainty
    '''
    from crpropa import Random, nucleusId, EeV, Vector3d, kpc, Candidate, ParticleState
    from useful_funcs import eqToGal
    import math

    R = Random(seed)  # CRPropa random number generator
    #Pid: 1000000000 + Z * 10000 + A * 10
    H_id = '-1000010010'
    He_id = '-1000020040'
    C_id = '-1000060120'
    N_id = '-1000070140'
    Fe_id = '-1000260520'

    all_events_lons, all_events_lats = [], []
    initial_lons, initial_lats = [], []
    NUM_OF_SIMS = 10
    for event in events:
        #if int(event[0]) != 3: continue
        particles = [- nucleusId(1,1), - nucleusId(4,2), - nucleusId(12,6), - nucleusId(52,26)]

        mean_energy = event[3] * EeV
        sigma_energy = 0.1 * mean_energy  # 10% energy uncertainty
        #sigma_energy = 0
        position = Vector3d(-8.5, 0, 0) * kpc

        lon0,lat0 = eqToGal(event[1], event[2])        #RETURN WHEN NO TEST
        #lon0,lat0 = event[1]*math.pi/180, event[2]*math.pi/180
        lat0 = math.pi/2 - lat0 #CrPropa uses colatitude, e.g. 90 - lat in degrees
        mean_dir = Vector3d()
        mean_dir.setRThetaPhi(1, lat0, lon0)
        sigma_dir = 0.002 # - 1 degree directional uncertainty
        initial_lons.append(lon0)
        initial_lats.append(lat0)

        lons_event, lats_event = [], []
        for pid in particles:
            lons, lats = [], []
            for i in range(NUM_OF_SIMS):
                energy = R.randNorm(mean_energy, sigma_energy)
                direction = R.randVectorAroundMean(mean_dir, sigma_dir)

                candidate = Candidate(ParticleState(pid, energy, position, direction))
                sim.run(candidate)

                res_direction = candidate.current.getDirection()
                #print(f"Event {event[0]}: {candidate}, direction:{res_direction}")
                lons.append(res_direction.getPhi())
                lats.append(res_direction.getTheta())
            if str(pid) == H_id:
                lons_event.append((lons, 'H'))
                lats_event.append((lats, 'H'))
            elif str(pid) == He_id:
                lons_event.append((lons, 'He'))
                lats_event.append((lats, 'He'))
            elif str(pid) == Fe_id:
                lons_event.append((lons, 'Fe'))
                lats_event.append((lats, 'Fe'))
            else:
                lons_event.append((lons, 'C'))
                lats_event.append((lats, 'C'))

        all_events_lats.extend(lats_event)
        all_events_lons.extend(lons_event)

    return initial_lats, initial_lons, all_events_lats, all_events_lons
