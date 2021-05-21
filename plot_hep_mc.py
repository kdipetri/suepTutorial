#!/usr/bin/env python

import pyhepmc_ng as hep
import matplotlib.pyplot as plt
import numpy as np
from particle import Particle
import math

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

        self.doTest = False 
        self.maxEvents = 100
        self.trackPtCut = 0.7
        self.maxEtaCut = 2.5 

        #
        # Used for plotting, add what you like
        # 
        self.higgsPt  = []
        self.higgsEta = []
        self.higgsPhi = []
        self.higgsM   = []

        self.nCharged   = [] 
        self.chargedPt  = []

        self.nTracks    = []
        self.trackPt     = []
        self.trackEta    = []
        self.trackPhi    = []
        self.trackM      = []

        self.isotropy   = []
        self.ht   = []

    def getDecayProducts(self,parent):
        # recursive function to get all status 1 decay products
        statusOne = []
    
        for particle in parent.children: 
          if particle.status == 1 : statusOne.append(particle)
          else : statusOne.extend( self.getDecayProducts(particle) ) 
    
        return statusOne

    def isCharged(self,particle):
        part = Particle.from_pdgid(particle.pid)
        if part.charge !=0 : return 1
        else : return 0
    
    def getIsotropy(self,particles,nSteps=1000):
        # https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/CandUtils/src/EventShapeVariables.cc
        # compute event isotropy
        deltaPhi = 2 * math.pi / nSteps  
        phi,eIn,eOut = 0,-1.,-1. 
    
        for i in range(0,nSteps):
            phi += deltaPhi
            tot = 0
            cosphi = math.cos(phi)
            sinphi = math.sin(phi)
            for particle in particles: 
                # sum over inner product of unit vectors and momenta
                particle_mom = particle.momentum
                particle_x = particle_mom.x
                particle_y = particle_mom.y
                tot += math.fabs( cosphi * particle_x + sinphi * particle_y )
            if eOut < 0. or tot < eOut : eOut = tot
            if eIn < 0. or tot > eIn : eIn = tot
    
        return 1.-(eIn - eOut)/eIn

    def processEvents(self):
        
        print( "Processing {}".format(self.infile) )

        # Reads the file
        with hep.open(self.infile) as f:
          # Just keeps looping
          while True :
        
            # Try to get an event
            evt = f.read()
            if not evt : break # we're at the end of the file
        
            # Stop if this is just a test or we're past max events
            if self.doTest and evt.event_number >= 5 :
                break
            if evt.event_number > self.maxEvents : 
                break
        
            # From here on, do things with the event!
            if self.doTest : 
                print("Event",evt.event_number)
         
            for particle in evt.particles :
            # This is what's in the "particle" class: http://hepmc.web.cern.ch/hepmc/classHepMC3_1_1GenParticle.html
        
              # Find the higgs
              if (particle.pid != 25 ) : continue 
              if (particle.status != 62 ) : continue
        
              # Get the higgs four vector
              evt_higgsPt   = particle.momentum.pt()
              evt_higgsEta  = particle.momentum.eta()
              evt_higgsPhi  = particle.momentum.phi()
              evt_higgsM    = particle.momentum.m()

              self.higgsPt .append(evt_higgsPt) 
              self.higgsEta.append(evt_higgsEta)
              self.higgsPhi.append(evt_higgsPhi)
              self.higgsM  .append(evt_higgsM)

              if self.doTest : 
                  print("higgs(pt,eta,phi,m)  = {:.1f},{:.2f},{:.2f},{:.1f}".format(evt_higgsPt,evt_higgsEta,evt_higgsPhi,evt_higgsM))

            
              # 
              # Study final state charged particles 
              # 
              decayProducts = self.getDecayProducts(particle)
        
              evt_ht = 0 # scalar sum track pt
              evt_ncharged = 0
              evt_ntracks = 0
              evt_tracks = []
              for child in decayProducts:

                  # only charged particles
                  if not self.isCharged(child) : continue 

                  evt_ncharged += 1
                  self.chargedPt.append(child.momentum.pt())
        
                  # charged particles in tracking acceptance
                  if child.momentum.pt() < self.trackPtCut : continue
                  if math.fabs(child.momentum.eta()) > self.maxEtaCut : continue

                  evt_ntracks += 1
                  evt_tracks.append(child)

                  self.trackPt .append(child.momentum.pt())
                  self.trackEta.append(child.momentum.eta())
                  self.trackPhi.append(child.momentum.phi())
                  self.trackM  .append(child.momentum.m())
                  evt_ht += child.momentum.pt()
              
              self.nCharged.append(evt_ncharged)
              self.nTracks .append(evt_ntracks)

              evt_isotropy = self.getIsotropy(evt_tracks)
              self.isotropy.append(evt_isotropy)

              self.ht.append(evt_ht)
        
              if self.doTest: 
                print("Number of charged particles:",evt_ncharged) 
                print("Isotropy: {:.2f}".format( evt_isotropy) )

        return

    def basicPlots(self,save=False): 
        # Plots basic histograms 
        fig, axs = plt.subplots(2,2,sharey=False,figsize=(10,8))
        #fig, axs = plt.subplots(2,2,sharey=False, tight_layout=True)
        
        # event display with track pT > 0.7 GeV
        nbins = 10
        axs[0,0].hist(self.nCharged  , bins=nbins )
        axs[0,0].set_xlabel("$n_{charged}$")
        axs[0,0].set_ylabel("events")
        
        # charged particle pT 
        axs[0,1].hist(self.chargedPt, bins=nbins*10)
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

# Testing the class
#suep = SUEP("mMed-400_mDark-1.0_temp-1.0_decay-generic.hepmc") 
#suep.maxEvents=100
#suep.doTest = False 
#suep.processEvents()
#suep.basicPlots()

