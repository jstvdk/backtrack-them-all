from useful_funcs import eqToGal
import typing
from run_sim import runSimulation
from visualizer import visualizeTotal, SimMap
from make_data import makeDF

def setupSimulation():
    '''
    This function setups current simulation
    '''
    # magnetic field setup
    B = JF12Field()
    seed = 691342
    B.randomStriated(seed)
    B.randomTurbulent(seed)
    #B = TF17Field()
    #B = PT11Field()

    # simulation setup
    sim = ModuleList()
    sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
    obs = Observer()
    obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
    sim.add(obs)
    #print(sim)

    return sim, obs

if __name__ == '__main__':
    import math
    from crpropa import *

    '''
    Getting data
    '''
    '''
    MAIN
    '''
    events = []
    with open('data/AugerApJS2022_highE.dat', 'r') as infile:
        for line in infile:
            if line.split()[0] == '#': continue
            temp_event = line.split()
            events.append((temp_event[0], float(temp_event[6]), float(temp_event[7]), float(temp_event[8])))

    '''
    tests
    '''
    '''
    events = []
    with open('data/test_data.csv', 'r') as infile:
        for line in infile:
            temp_event = line.split(',')
            events.append((temp_event[0], float(temp_event[1]), float(temp_event[2]), float(temp_event[3])))
    '''
    #events = []
    #with open('data/nesessary_objects.csv', 'r') as infile:
    #    for line in infile:
    #        temp_event = line.split(',')
    #        events.append((temp_event[0], float(temp_event[1]), float(temp_event[2]), float(temp_event[3])))
    '''
    Setupping simulation
    '''
    sim, obs = setupSimulation()

    '''
    Running simulation
    '''
    initial_lats, initial_lons, all_events_lats, all_events_lons = runSimulation(sim, obs, events, seed=42)

    '''
    GATHERING DATA
    '''

    total_results = makeDF(all_events_lats, all_events_lons)

    '''
    Visualizing results achieved
    '''
    map = SimMap(total_results, initial_lats, initial_lons, particles=['Fe'])
    map.setSaveName('results/full_ehecr_map/Fe_backtrack.png')
    map.setTitle("Simulated CRs with E > 20 Eev for Z = 26")
    map.setSourcesFlags({'mags': True, 'sbgs': True, 'clusts': True})
    map.plotMap(transform=True)
    #visualizeTotal(total_results, initial_lats, initial_lons)


    '''
    Saving results
    '''
    #total_results.to_csv('auger31_42.csv')
