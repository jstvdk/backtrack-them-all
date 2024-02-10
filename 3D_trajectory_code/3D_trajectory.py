from crpropa import *
from useful_funcs import eqToGal
import math


class MyTrajectoryOutput(Module):
    """
    Custom trajectory output: i, x, y, z
    where i is a running cosmic ray number
    and x,y,z are the Galactocentric coordinates in [kpc].
    """
    def __init__(self, fname):
        Module.__init__(self)
        self.fout = open(fname, 'w')
        self.fout.write('#i\tX\tY\tZ\n')
        self.i = 0
    def process(self, c):
        v = c.current.getPosition()
        x = v.x / kpc
        y = v.y / kpc
        z = v.z / kpc
        self.fout.write('%i\t%.3f\t%.3f\t%.3f\n'%(self.i, x, y, z))
        if not(c.isActive()):
            self.i += 1
    def close(self):
        self.fout.close()


if __name__ == '__main__':
    # magnetic field setup
    seed = 691342
    R = Random(seed)
    B = JF12Field()
    B.randomStriated(seed)
    B.randomTurbulent(seed)

    # simulation setup
    sim = ModuleList()
    sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
    sim.add(SphericalBoundary(Vector3d(0), 20 * kpc))
    output = MyTrajectoryOutput('galactic_trajectories_with_uncert_with_random_4types.txt')
    sim.add(output)
    NUM_OF_SIMS = 1

    '''
    Creating events list
    '''
    events = []
    with open('data/AugerApJS2022_highE.dat', 'r') as infile:
        for line in infile:
            if line.split()[0] == '#': continue
            temp_event = line.split()
            events.append((temp_event[0], float(temp_event[6]), float(temp_event[7]), float(temp_event[8])))

    '''
    Sim for 4 particles for 1 event(third one)
    '''
    #particles = [- nucleusId(1,1), - nucleusId(4,2), - nucleusId(12,6), - nucleusId(52,26)]
    particles = [- nucleusId(1,1), - nucleusId(4,2), - nucleusId(12,6), - nucleusId(52,26)]
    for event in events:
        if int(event[0]) != 3: continue

        meanEnergy = event[3] * EeV
        sigmaEnergy = 0.1 * meanEnergy  # 10% energy uncertainty
        #sigmaEnergy = 0
        position = Vector3d(-8.5, 0, 0) * kpc

        lon0,lat0 = eqToGal(event[1], event[2])        #RETURN WHEN NO TEST
        lat0 = math.pi/2 - lat0 #CrPropa uses colatitude, e.g. 90 - lat in degrees
        meanDir = Vector3d()
        meanDir.setRThetaPhi(1, lat0, lon0)
        #sigmaDir = 0.0
        sigmaDir = 0.002 #1 degree directional uncertainty

        for pid in particles:
            for i in range(NUM_OF_SIMS):
                energy = R.randNorm(meanEnergy, sigmaEnergy)
                direction = R.randVectorAroundMean(meanDir, sigmaDir)

                candidate = Candidate(ParticleState(pid, energy, position, direction))
                sim.run(candidate)

    output.close()
