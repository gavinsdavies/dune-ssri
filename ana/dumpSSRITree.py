#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
import glob

muon_mass = 105.6583755
extra_trk_length = 24.35 # g/cm^2 from Mike K.
nmax = -1
maxfiles = 9999

def loop( events, dspt, tgeo, tout ):

    offset = [ 0., 5.5, 411. ]
    fvLo = [ -300., -100., 50. ]
    fvHi = [ 300., 100., 450. ]
    collarLo = [ -320., -120., 30. ]
    collarHi = [ 320., 120., 470. ]

    nplanes_cut = 8

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    print "Starting loop over %d entries" % N
    ient = 0
    for evt in dspt:

        t_Ev[0] = evt.StdHepP4[3] # neutrino is first 4-vector, so 3 is always Ev

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)

        if ient > nmax: 
          break

        events.GetEntry(ient)

        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            t_ievt[0] = ient
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            # In SSRI
            t_muonDeath[0]=0.0; t_muonDeath[1]=0.0; t_muonDeath[2]=0.0;
            t_muonBirth[0]=0.0; t_muonBirth[1]=0.0; t_muonBirth[2]=0.0;

            t_lepTMS[0]=0.0; t_lepTMS[1]=0.0; t_lepTMS[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_TMSKE[0] = -1.
            t_lepE[0] = 0.
            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonExitKE[0] = 0.0
            t_muonReco[0] = -1;
            t_muScintLen[0] = 0.0;
            t_muLArLen[0] = 0.0;
            t_muScintLenOld[0] = 0.0;
            t_muScintEnergy[0] = 0.0;

            xpt.clear()
            zpt.clear()
            deposit.clear()

            # now ID numucc
            t_Reac=vertex.GetReaction()

            # set the vertex location for output
            for i in range(3): t_vtx[i] = vertex.GetPosition()[i] / 10. - offset[i] # cm

            ileptraj = -1
            nfsp = 0
            # get the lepton kinematics from the edepsim file
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.GetMomentum()[3]
                # Get the momentum (px^2+py^2+pz^2)
                p = (particle.GetMomentum()[0]**2 + particle.GetMomentum()[1]**2 + particle.GetMomentum()[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.GetPDGCode()
                t_fsPx[nfsp] = particle.GetMomentum()[0]
                t_fsPy[nfsp] = particle.GetMomentum()[1]
                t_fsPz[nfsp] = particle.GetMomentum()[2]
                t_fsE[nfsp] = e
                nfsp += 1
                pdg = particle.GetPDGCode()
                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.GetTrackId()
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.GetMomentum()[i]
                    t_lepKE[0] = e - m
                    t_lepE[0] = e

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp

            if abs(t_lepPdg[0]) != 13: continue

            # If there is a muon, determine how to reconstruct its momentum and charge
            exit = False
            inTMS = False

            leptraj = event.Trajectories[ileptraj]

            pPrev = None
            pointval = 0
            for p in leptraj.Points:
                pt = p.GetPosition()
                node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                volName = node.GetName()
                active = False

                pPos = ROOT.TVector3( pt.X()/10. - offset[0], pt.Y()/10. - offset[1], pt.Z()/10. - offset[2] )
                #print "point ",pointval,pPos.z(),volName

                # The z-information is contained in thinlayerposition0, thinlayervol

                if "LAr" in volName or "PixelPlane" in volName or "sPlane" in volName: # in active volume, update exit points
                    t_muonExitPt[0] = pt.X() / 10. - offset[0]
                    t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_muonExitMom[0] = p.GetMomentum().x()
                    t_muonExitMom[1] = p.GetMomentum().y()
                    t_muonExitMom[2] = p.GetMomentum().z()
                    
                    t_muonExitKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                else:
                    if not exit:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        exit = True

                    # Check if it is in the TMS
                    if ("TMS" in volName or "modulelayer" in volName) and not inTMS:
                        t_TMSKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                        inTMS = True
                pPrev = pPos
                pointval += 1

            endpt = leptraj.Points[-1].GetPosition()

            node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

            t_lepDeath[0] = endpt.X()/10. - offset[0]
            t_lepDeath[1] = endpt.Y()/10. - offset[1]
            t_lepDeath[2] = endpt.Z()/10. - offset[2]

            endVolName = node.GetName()
            if "TPCActive" in endVolName: t_muonReco[0] = 1 # contained
            # Updated name in new geom
            elif "TMS" in endVolName: t_muonReco[0] = 2 # Scintillator stopper
            else: t_muonReco[0] = 0 # endpoint not in active material

            # look for muon hits in the ArgonCube
            arhits = []
            for key in event.SegmentDetectors:
                # Updated name in new geom
                if key.first == "TPCActive_shape":
                    arhits += key.second
                    
            ar_muon_hits = []
            for idx, hit in enumerate(arhits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                #if traj.GetParentId() == -1 and abs(traj.GetPDGCode()) == 13:
                ar_muon_hits.append(hit)

            if len(ar_muon_hits) < 3: continue

            ar_trk_length_gcm2 = 0.
            for idx, hit in enumerate(ar_muon_hits):
                hStart = ROOT.TVector3( hit.GetStart()[0]/10.-offset[0], hit.GetStart()[1]/10.-offset[1], hit.GetStart()[2]/10.-offset[2] )
                xpt.push_back(hStart.x())
                zpt.push_back(hStart.z())
                
            # just use vertex and LAr endpoint for first pass
            hArStart=ROOT.TVector3(ar_muon_hits[0].GetStart()[0]/10.-offset[0], ar_muon_hits[0].GetStart()[1]/10.-offset[1], ar_muon_hits[0].GetStart()[2]/10.-offset[2])
            hArEnd=ROOT.TVector3(ar_muon_hits[-1].GetStart()[0]/10.-offset[0], ar_muon_hits[-1].GetStart()[1]/10.-offset[1], ar_muon_hits[-1].GetStart()[2]/10.-offset[2])

            ar_dh = (hArEnd-hArStart).Mag()
            ar_dz = hArEnd.z()-hArStart.z()
                
            ar_trk_length_gcm2 = 1.4 * ar_dh

            #-------------------------------------------------------
            # look for muon hits in the scintillator
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "volTMS":
                    hits += key.second

            muon_hits = []
            for idx, hit in enumerate(hits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                #if traj.GetParentId() == -1 and abs(traj.GetPDGCode()) == 13:
                muon_hits.append(hit)

            if len(muon_hits) < 2: continue

            hMuonStart = muon_hits[0].GetStart()
            t_muonBirth[0] = hMuonStart[0]/10.-offset[0]
            t_muonBirth[1] = hMuonStart[1]/10.-offset[1]
            t_muonBirth[2] = hMuonStart[2]/10.-offset[2]

            trk_length_gcm2 = 0.
            de = 0.

            hPrev = None
            hFinal = None
            hPrevTime = 0
            total_planes = 0
            for idx, hit in enumerate(muon_hits):
                de += hit.GetEnergyDeposit()

                hStart = ROOT.TVector3( hit.GetStart()[0]/10.-offset[0], hit.GetStart()[1]/10.-offset[1], hit.GetStart()[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.GetStop()[0]/10.-offset[0], hit.GetStop()[1]/10.-offset[1], hit.GetStop()[2]/10.-offset[2] )
                Time = hit.GetStart()[3]

                #tempnode = tgeo.FindNode( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2])
                #nav = tgeo.GetCurrentNavigator()
                #topnode_name = nav.GetCurrentNode().GetName()
                # Z enumerated 
                #while "modulelayervol_PV" not in topnode_name:
                  #nav.CdUp()
                  #topnode_name = nav.GetCurrentNode().GetName()

                #print topnode_name
                #layerno = nav.GetCurrentNode().GetNumber()
                #print layerno
                #hStart.Print()

                #nextnode = nav.FindNextBoundary();
                #print "Next boundary: ",nextnode.GetName(), nextnode.GetVolume().GetName()

                xpt.push_back(hStart.x())
                zpt.push_back(hStart.z())
                deposit.push_back(hit.GetEnergyDeposit())
                continue

                # this isn't the first hit, so we can start to build a track
                if hPrev is not None:

                    # thin layer is 3cm air + 1cm scint + 1.5cm steel = 5.5cm pitch
                    # -> means 3cm+1.5cm = 4.5cm between hits
                    # Runs from 733 to 949

                    # thick layer is 3cm air + 1cm scint + 4cm steel = 8cm pitch
                    # -> means 3cm+4cm = 7cm between hits

                    # There is also the transition layer
                    # 1.5cm steel + 4cm steel stuck together, with air on either side (each 1.5cm) 
                    # -> 1.5+4cm+1.5+1.5 = 8.5cm between hits
                    # find 8.5cm between hits

                    # if "gap" isn't a multiple of that, then something has gone off the rails and we want to stop
                    gap = hStart.z() - hPrev.z()
                    # correct for cosine of the incident angle
                    dh = (hStart-hPrev).Mag()
                    dz = hStart.z()-hPrev.z()

                    #Found change at 730.8-736.3
                    #Previous distance: 730.8
                    #New distance: 5.5

                    #Found change at 939.8-949.3
                    #Previous distance: 5.5
                    #New distance: 9.5

                    #Found change at 949.3-957.3
                    #Previous distance: 9.5
                    #New distance: 8

                    
                    # sometimes there are multiple hits within a scintillator plane, and we want to skip them
                    # The smallest gap between hits is 3cm air +1.5cm steel -> 4.5cm

                    # Inspect the geometry, find if this position is the same node
                    same = tgeo.IsSameLocation( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2])
                    tempnode = tgeo.FindNode( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2])
                    if same == True:
                      continue
                    if gap < 4.5:
                      continue

                    timediff = Time-hPrevTime
                    # Make a cut on the time difference to avoid delayed hits
                    # Only when the xhit is not in the dead region
                    if timediff > 2 and (abs(hStart.x()) > 182 or (abs(hStart.x()) < 170 and abs(hStart.x()) > 8)):
                      continue

                    # "thin layer" region up until 949.2925 (scintillator here), so gap should be 5.5
                    if hStart.z() < 949.:
                        nplanes = int(gap/5.5)
                        # if nplanes isn't an integer, then we've reached the end of the track
                        # and there is some other hit that is super far away typically that is messing up the length

                        total_planes += nplanes
                        # Sometimes the track goes out of scintillator, still propagates through the iron and air, then comes back into the scintillator
                        # This leads to a larger gap, but it only goes through the iron, so modify track length
                        # The dh/dz modification (costheta) will also be borked (have no info on bend), so assume straight line
                        if gap > 5.5 and nplanes > 1:
                          trk_length_gcm2 += (1.5*7.85) * nplanes * dh/dz
                        else:
                          # we must have gone through 1cm scint + 1.5cm steel + 3cm air, divide by cos(theta)
                          trk_length_gcm2 += (1.05 + 1.5*7.85) * nplanes * dh/dz

                        '''
                        if abs(nplanes - gap/5.5) > 1.E-3: 
                          print nplanes
                          print abs(nplanes-gap/5.5), gap, ient
                          hPrev.Print()
                          hStart.Print()
                          break
                        if nplanes > nplanes_cut: # find an intolerable gap?
                            break
                        '''


                    # The transition layer gap
                    # Between 939.8+1 to 949.3
                    #elif (hStart.z() - 949.3) < 0.1: # the transition module, the gap is 9.5
                    elif (hStart.z() < 949.3+1.5): # Scintillator strip goes from 949.3 to 949.3+1 (add in 1.5 for good measure)
                        if abs(gap-9.5) < 1.E-3: # Gap is 9.5cm
                            trk_length_gcm2 += (1.05 + (4.+1.5)*7.85) * dh/dz # Gone through 1 cm scintillator, 4cm steel, 1.5cm steel, 3cm air
                            total_planes += 1 # Add in the plane

                    else: #"thick layer region, gap should be 8cm
                        nplanes = int(gap/8.)
                        '''
                        if abs(nplanes - gap/8.) > 1.E-3: # the track is gone to shit
                            print nplanes
                            print abs(nplanes-gap/8.), gap, ient
                            hPrev.Print()
                            hStart.Print()
                            break
                        if nplanes > nplanes_cut:
                            break
                        '''
                        total_planes += nplanes

                        if gap > 8 and nplanes > 1:
                          trk_length_gcm2 += (4.*7.85) * nplanes * dh/dz
                        else:
                          trk_length_gcm2 += (1.05 + 4.*7.85) * nplanes * dh/dz
                hPrev = hStart # update it
                hPrevTime = Time
                hFinal = hStop

            hit = muon_hits[-1]
            hFinal = ROOT.TVector3( hit.GetStart()[0]/10.-offset[0], hit.GetStart()[1]/10.-offset[1], hit.GetStart()[2]/10.-offset[2] )
            t_muonDeath[0] = hFinal.x()
            t_muonDeath[1] = hFinal.y()
            t_muonDeath[2] = hFinal.z()

            #t_muScintLen[0] = trk_length_gcm2 + extra_trk_length
            # Take the directionality into account when adding the track length
            t_muScintLen[0] = trk_length_gcm2 + extra_trk_length
            t_muLArLen[0] = ar_trk_length_gcm2 
            t_muScintEnergy[0] = de
            
            tout.Fill()
        ient += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    #parser.add_option('--outfile', help='Output file name', default="out.root")
    #parser.add_option('--topdir', help='Input file top directory', default="")
    #parser.add_option('--nmax', help='Maximum number of events', default="99999999999999")
    #parser.add_option('--maxfiles', help='Maximum number of events', default="99999")
    #parser.add_option('--first_run', type=int, help='First run number', default=0)
    #parser.add_option('--last_run', type=int, help='Last run number', default=0)
    #parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    #parser.add_option('--grid', action='store_true', help='grid mode', default=False)
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--nmax', help='Maximum number of events', default="99999999999999")
    parser.add_option('--infile', help='EDepSim input file')

    (args, dummy) = parser.parse_args()

    nmax = int(args.nmax)
    input_file = args.infile
    print input_file
    if not input_file:
      print "Need input file"
      exit(-1)

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_Ev = array('f', [0.])
    tout.Branch('Ev',t_Ev,'Ev/F')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')
    t_muonDeath = array('f',3*[0.0])
    tout.Branch('muonDeath',t_muonDeath,'muonDeath[3]/F')
    t_muonBirth = array('f',3*[0.0])
    tout.Branch('muonBirth',t_muonBirth,'muonBirth[3]/F')
    t_lepTMS = array('f',3*[0.0])
    tout.Branch('lepTMS',t_lepTMS,'leTMS[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_TMSKE = array('f',[0])
    tout.Branch('TMSKE',t_TMSKE,'TMSKE/F')
    t_lepE = array('f',[0])
    tout.Branch('lepE',t_lepE,'lepE/F')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muScintLen = array('f',[0])
    tout.Branch('muScintLen',t_muScintLen,'muScintLen/F')
    t_muLArLen = array('f',[0])
    tout.Branch('muLArLen',t_muLArLen,'muLArLen/F')
    t_muonExitKE = array('f',[0])
    tout.Branch('muonExitKE',t_muonExitKE,'muonExitKE/F')
    t_muScintLenOld = array('f',[0])
    tout.Branch('muScintLenOld',t_muScintLenOld,'muScintLenOld/F')
    t_muScintEnergy = array('f',[0])
    tout.Branch('muScintEnergy',t_muScintEnergy,'muScintEnergy/F')
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')

    t_Reac = "";
    tout.Branch('reac', t_Reac, 'reac/C')

    xpt = ROOT.std.vector('float')()
    zpt = ROOT.std.vector('float')()
    deposit = ROOT.std.vector('float')()
    tout.Branch('xpt', xpt)
    tout.Branch('zpt', zpt)
    tout.Branch('deposit', deposit)
    
    loaded = False
    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    #neutrino = "neutrino"
    #if args.rhc:
        #neutrino = "antineutrino"

    #print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    #for run in range( args.first_run, args.last_run+1 ):
        #if args.grid:
            #fname = "%s/edep.%d.root" % (args.topdir,run)
        #else:
            #fname = "%s/%s.%d.edepsim.root" % (args.topdir, neutrino, run)
        #print fname
#
        #fname = fname.replace("/pnfs","root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr")
        # see if it is an OK file
        #if not os.access( fname, os.R_OK ):
            #print "Can't access file: %s" % fname
            #continue
        #tf = ROOT.TFile.Open( fname )
        #if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            #print "File is crap: %s" % fname
            #continue
#
        #if not loaded:
            #loaded = True
            #tf.MakeProject("EDepSimEvents","*","RECREATE++")
##
        # add it to the tchain
        #events.Add( fname )
        #dspt.Add( fname )

        #if tgeo is None: # first OK file, get geometry
            #tgeo = tf.Get("EDepSimGeometry")
        #tf.Close() # done with this one

    #nfiles = 0
    #for run in glob.glob(args.topdir+"/edep*.root"):
      #print nfiles
      #print maxfiles
      #if nfiles > maxfiles: 
        #break

      #fname = run
    fname = input_file
    print "Adding "+fname+" to TChain..."
    tgeo = None

    tf = ROOT.TFile.Open( fname )
    tf.MakeProject("EDepSimEvents","*","RECREATE++")

    # add it to the tchain
    events.Add( fname )
    dspt.Add( fname )

    #if tgeo is None: # first OK file, get geometry
    tgeo = tf.Get("EDepSimGeometry")
    tf.Close() # done with this one

    #nfiles += 1

    #print "Running on ",nfiles," root files..."
    loop( events, dspt, tgeo, tout )

    fout.cd()
    tout.Write()

