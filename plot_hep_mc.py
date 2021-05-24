#!/usr/bin/env python
# based Katherine Pachal & Jess Nelson's skeleton 
# https://gitlab.cern.ch/snowmass-track-trigger/scripts_lhe_hepmc
# event displays from Chris Papageorgakis

from particle import Particle
import pyhepmc_ng as hep
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

def getStatusOne(particles):
    # get all status 1 decay products
    statusOne = []

    for particle in particles: 
      if particle.status != 1 : continue
      statusOne.append(particle)

    return statusOne

def hasScalarParent(particle):
    # returns True if a particle descends from scalar 
    scalar_parent=False

    parents=particle.parents
    if not parents : return scalar_parent
    else:
        for parent in particle.parents:
            if parent.pid == 999999 :return True# dark meson
            else : return hasScalarParent(parent)

    return statusOne

def isCharged(particle):
    part = Particle.from_pdgid(particle.pid)
    if part.charge !=0 : return 1
    else : return 0

def getIsotropy(particles,nSteps=100):
    # compute event isotropy, based on.... 
    # https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/CandUtils/src/EventShapeVariables.cc

    def getSum(phi,particle):
        cosphi, sinphi = math.cos(phi), math.sin(phi)
        return math.fabs( cosphi * particle + sinphi * particle )

    phis = np.linspace(0,2*math.pi,nSteps)
    cosphis = np.cos(phis)
    sinphis = np.sin(phis)
    sums = np.zeros(nSteps)
    for particle in particles:
        sums += np.abs( cosphis * particle.momentum.x + sinphis * particle.momentum.y)

    eOut = min(sums)
    eIn = max(sums)

    return 1.-(eIn - eOut)/eIn


class SUEP():

    def __init__(self, infile):
        #
        # Get mediator mass, dark meson mass, temp, etc from filename
        # 
        self.infile = infile
        self.mMed  = infile.split("_")[0].split("-")[1] 
        self.mDark = infile.split("_")[1].split("-")[1] 
        self.temp  = infile.split("_")[2].split("-")[1] 
        self.decay = infile.split("_")[3].split("-")[1] 

        # Configurable options for running 
        self.verbose = False
        self.maxEvents = 100
        self.trackPtCut = 0.7
        self.maxEtaCut = 2.5 

        #
        # Used for plotting, add what you like
        # 
        self.scalarPt  = []
        self.scalarEta = []
        self.scalarPhi = []
        self.scalarM   = []

        self.nCharged   = [] 
        self.chargedPt  = []

        self.nTracks    = []
        self.trackPt     = []
        self.trackEta    = []
        self.trackPhi    = []

        self.isotropy   = []
        self.ht   = []


    def getScalarInfo(self,particles):
        # Find the mediator
        for particle in particles :
            
            # Find the scalar
            if (particle.pid != 25 ) : continue 
            if (particle.status != 62 ) : continue

            # Get the scalar four vector

            self.scalarPt .append(particle.momentum.pt()) 
            self.scalarEta.append(particle.momentum.eta())
            self.scalarPhi.append(particle.momentum.phi())
            self.scalarM  .append(particle.momentum.m())
            
            break

        return

    def saveTrackInfo(self,particles):

        evt_ht = 0 # scalar sum track pt
        evt_ncharged = 0 
        evt_ntracks = 0
        evt_tracks = []
        evt_trackEta = []
        evt_trackPhi = []
        evt_trackPt  = []

        for particle in getStatusOne(particles):

            # only charged particles
            if not isCharged(particle) : continue 

            evt_ncharged += 1
            self.chargedPt.append(particle.momentum.pt())
        
            # charged particles in tracking acceptance
            if particle.momentum.pt() < self.trackPtCut : continue
            if math.fabs(particle.momentum.eta()) > self.maxEtaCut : continue

            evt_ntracks += 1
            evt_tracks.append(particle)

            evt_ht += particle.momentum.pt()
            evt_trackPt .append(particle.momentum.pt())
            evt_trackEta.append(particle.momentum.eta())
            evt_trackPhi.append(particle.momentum.phi())

        # for event displays
        self.trackPt .append(evt_trackPt)
        self.trackEta.append(evt_trackEta)
        self.trackPhi.append(evt_trackPhi)
        
        # Event level quantities 
        self.nCharged.append(evt_ncharged)
        self.nTracks .append(evt_ntracks)
        self.ht.append(evt_ht)
        self.isotropy.append( getIsotropy(evt_tracks) )

        return 

    def processEvents(self):
        
        print( "Processing {}".format(self.infile) )

        # Reads the file
        with hep.open(self.infile) as f:
          # Just keeps looping
          while True :
        
            # Try to get an event
            evt = f.read()
            if not evt : break # we're at the end of the file
        
            # Stop if we're past max events
            if evt.event_number >= self.maxEvents : 
                break
        
            # Study the mediator and final state charged particles 
            self.getScalarInfo(evt.particles)
            # End looping because we found the scalar
            self.saveTrackInfo(evt.particles)

            if self.verbose : 
                print("Event",evt.event_number)

                # Mediator info
                print("scalar: PtEtaPhiM[{:.2f}, {:.2f}, {:.2f}, {:.2f}]".format(self.scalarPt[-1],self.scalarEta[-1],self.scalarPhi[-1],self.scalarM[-1]) )

                # Status=1 particles
                print("Some final state charged particles: ")
                particles = getStatusOne(evt.particles)
                i=0
                for part in particles : 
                  if not hasScalarParent(part): continue
                  if part.momentum.pt() < 0.7: continue 
                  if not isCharged(part): continue
                  i+=1
                  trk=part.momentum
                  print("pdgID {} : PtEtaPhiM[{:.2f}, {:.2f}, {:.2f}, {:.2f}]".format(part.pid,trk.pt(),trk.eta(),trk.phi(),trk.m()))
                  if i>15: break
                # Event level info
                print("Number of status=1 particles:",len(particles)) 
                print("Number of charged particles:",self.nCharged[-1]) 
                print("Number of truth tracks:",self.nTracks[-1]) 
                print("Isotropy: {:.2f}".format(self.isotropy[-1] ))
                print("HT: {:.1f}".format(self.ht[-1] ))
                

        return

    def doTest(self):
        self.verbose=True
        self.maxEvents=1
        self.processEvents()
        self.verbose=False
        self.maxEvents=100
        return
    
    def basicPlots(self,save=False): 
        # Plots basic histograms 
        fig, axs = plt.subplots(2,2,sharey=False,figsize=(10,8))
        
        # number of status=1 charged particles 
        nbins = 10
        axs[0,0].hist(self.nCharged  , bins=nbins )
        axs[0,0].set_xlabel("$n_{charged}$")
        axs[0,0].set_ylabel("events")
        
        # charged particle pT 
        axs[0,1].hist(self.chargedPt, bins=np.linspace(0,5,50))
        axs[0,1].set_xlim(0,5)
        axs[0,1].set_xlabel("$p_{T}^{charged}$")
        axs[0,1].set_ylabel("charged particles")
        
        # n truth tracks pT > 0.7 GeV & |eta|<2.5
        axs[1,0].hist(self.nTracks  , bins=nbins)
        axs[1,0].set_xlabel("$n_{tracks}$")
        axs[1,0].set_ylabel("events")
        
        # isotropy with truth tracks 
        axs[1,1].hist(self.isotropy , bins=nbins)
        axs[1,1].set_xlabel("isotropy")
        axs[1,1].set_ylabel("events")
        
        #plt.show()
        if (save) : plt.savefig("basicPlot_{}.png".format(self.infile.strip(".hepmc")))
        else : plt.show()
        return

    def eventDisplay(self,event=0,nbins=25,save=False): 
        # 
        # Makes SUEP event display 

        # quantites for plotting
        part_pt = np.array([],dtype=float)
        part_eta= np.array([],dtype=float)
        part_phi= np.array([],dtype=float)
        ht=0
        # masks for plotting
        fromScalar = np.array([],dtype=bool)
        fromIsr    = np.array([],dtype=bool)
        is_mu      = np.array([],dtype=bool)
        is_e       = np.array([],dtype=bool)
        is_gam     = np.array([],dtype=bool)
        is_had     = np.array([],dtype=bool)

        # Read the file
        with hep.open(self.infile) as f:
          while True :
            evt = f.read()
            if not evt : break # we're at the end of the file
            if evt.event_number!=event: continue
            
            particles = getStatusOne(evt.particles)
            for particle in particles:

                # get track info
                mom = particle.momentum
                #if mom.pt() < 0.1 : continue
                if mom.pt() < self.trackPtCut : continue
                #if abs(mom.eta()) > self.maxEtaCut : continue
                part_pt = np.append(part_pt , mom.pt())
                part_phi= np.append(part_phi, mom.phi())
                part_eta= np.append(part_eta, mom.eta())
                ht+=mom.pt()

                # handle masks
                parent = hasScalarParent(particle)
                fromScalar=np.append(fromScalar,parent)
                fromIsr   =np.append(fromIsr,not parent)
                is_e  =np.append(is_e  , abs(particle.pid) == 11)
                is_mu =np.append(is_mu , abs(particle.pid) == 13)
                is_gam=np.append(is_gam, abs(particle.pid) == 22)
                is_had=np.append(is_had, abs(particle.pid) > 100)

            break # done

        # 2D scatter plot 
        fig, ax = plt.subplots(figsize =(8, 6))
        def scale(vec_pt):
            s = [50.*pt for pt in vec_pt ]
            #s = [2500.*pt/ht for pt in vec_pt ]
            return s
        def draw(mask,mark='o',col='xkcd:blue'):
            ax.scatter(part_eta[mask], part_phi[mask], s=scale(part_pt[mask]),marker=mark,c=col)
            return
        
        draw(fromScalar & is_e  ,mark='o',col='xkcd:blue')#
        draw(fromScalar & is_mu ,mark='v',col='xkcd:blue')#
        draw(fromScalar & is_gam,mark='s',col='xkcd:blue')#
        draw(fromScalar & is_had,mark='*',col='xkcd:blue')#
        draw(fromIsr & is_e  ,mark='o',col='xkcd:magenta')#
        draw(fromIsr & is_mu ,mark='v',col='xkcd:magenta')#
        draw(fromIsr & is_gam,mark='s',col='xkcd:magenta')#
        draw(fromIsr & is_had,mark='*',col='xkcd:magenta')#
    
        ax.set_xlabel("$\eta_{particle}$")
        ax.set_ylabel("$\phi_{particle}$")
        ax.set_xlim(-4.0,4.0)
        ax.set_ylim(-3.5,3.5)

        # Legend 1 is particle type
        line1 = ax.scatter([-100], [-100], label='$e$',marker='o', c='xkcd:black')
        line2 = ax.scatter([-100], [-100], label='$\mu$', marker='v', c='xkcd:black')
        line3 = ax.scatter([-100], [-100], label='$\gamma$', marker='s', c='xkcd:black')
        line4 = ax.scatter([-100], [-100], label='hadron', marker='*', c='xkcd:black')
        first_legend = plt.legend(handles=[line1, line2, line3, line4],
                                  loc='upper right', fontsize=12)
        ax.add_artist(first_legend)
    
        # Legend 2 is particle origin
        blue_patch = mpatches.Patch(color='xkcd:blue', label='from scalar')
        red_patch  = mpatches.Patch(color='xkcd:magenta', label='not from scalar')
        plt.legend(handles=[blue_patch, red_patch],loc='upper left')
        
        if (save) : plt.savefig("eventdisplay_{}_{}.png".format(event,self.infile.strip(".hepmc")))
        else : plt.show()
        return

