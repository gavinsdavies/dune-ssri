#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array

muon_mass = 105.6583755
extra_trk_length = 24.35 # g/cm^2 from Mike K.

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

            t_lepRMMS[0]=0.0; t_lepRMMS[1]=0.0; t_lepRMMS[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_rmmsKE[0] = -1.
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

            #if ient != 5819: break
            # now ID numucc
            reaction=vertex.Reaction
            #print reaction

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm

            # fiducial vertex pre-cut
            #if abs(t_vtx[0]) > 310. or abs(t_vtx[1]) > 110. or t_vtx[2] < 40. or t_vtx[2] > 360.:
            #    continue

            ileptraj = -1
            nfsp = 0
            # get the lepton kinematics from the edepsim file
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                nfsp += 1
                pdg = particle.PDGCode
                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    #t_lepKE[0] = (particle.Momentum.Mag2() + muon_mass*2)**0.5 - muon_mass
                    #print "Lepton KE: ", t_lepKE[0]
                    t_lepKE[0] = e - m
                    t_lepE[0] = e
                    #print "Lepton KE: ", t_lepKE[0]

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp


            if abs(t_lepPdg[0]) != 13: continue

            # If there is a muon, determine how to reconstruct its momentum and charge
            exit = False
            inrmms = False

            leptraj = event.Trajectories[ileptraj]

            pPrev = None
            for p in leptraj.Points:
                pt = p.Position
                node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                volName = node.GetName()
                active = False

                pPos = ROOT.TVector3( pt.X()/10. - offset[0], pt.Y()/10. - offset[1], pt.Z()/10. - offset[2] )

                #if pPrev is not None:
                #    print "z %1.2f dz %1.2f vol %s" % (pPos.z(), pPos.z() - pPrev.z(), volName)
            
                if "LAr" in volName or "PixelPlane" in volName or "sPlane" in volName: # in active volume, update exit points
                    t_muonExitPt[0] = pt.X() / 10. - offset[0]
                    t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_muonExitMom[0] = p.Momentum.x()
                    t_muonExitMom[1] = p.Momentum.y()
                    t_muonExitMom[2] = p.Momentum.z()
                    
                    t_muonExitKE[0] = (p.Momentum.Mag2() + muon_mass*2)**0.5 - muon_mass
                    #print "Muon KE ", t_muonExitKE[0]
                else:
                    if not exit:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        exit = True

                    # Check if it is in the RMMS
                    if ("RMMS" in volName or "modulelayer" in volName) and not inrmms:
                        t_rmmsKE[0] = (p.Momentum.Mag2() + muon_mass*2)**0.5 - muon_mass
                        inrmms = True
                pPrev = pPos

            endpt = leptraj.Points[-1].Position

            node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

            t_lepDeath[0] = endpt.X()/10. - offset[0]
            t_lepDeath[1] = endpt.Y()/10. - offset[1]
            t_lepDeath[2] = endpt.Z()/10. - offset[2]

            endVolName = node.GetName()
            if "LArActive" in endVolName: t_muonReco[0] = 1 # contained
            elif "RMMS" in endVolName: t_muonReco[0] = 2 # Scintillator stopper
            else: t_muonReco[0] = 0 # endpoint not in active material

            # look for muon hits in the ArgonCube
            arhits = []
            for key in event.SegmentDetectors:
                if key.first == "ArgonCube":
                    arhits += key.second
                    
            ar_muon_hits = []
            for idx, hit in enumerate(arhits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                if traj.ParentId == -1 and abs(traj.PDGCode) == 13:
                    ar_muon_hits.append(hit)

            if len(ar_muon_hits) < 3: continue

            ar_trk_length_gcm2 = 0.
            #print ""
            for idx, hit in enumerate(ar_muon_hits):
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                #print "Start X: ", hStart.x()
                #print "Start Z: ", hStart.z()
                xpt.push_back(hStart.x())
                zpt.push_back(hStart.z())
                
            
            # just use vertex and LAr endpoint for first pass
            #hArStart = ROOT.TVector3 (t_vtx[0], t_vtx[1], t_vtx[2])
            #hArEnd = ROOT.TVector3( t_muonExitPt[0], t_muonExitPt[1], t_muonExitPt[2] )
            hArStart=ROOT.TVector3(ar_muon_hits[0].Start[0]/10.-offset[0], ar_muon_hits[0].Start[1]/10.-offset[1], ar_muon_hits[0].Start[2]/10.-offset[2])
            hArEnd=ROOT.TVector3(ar_muon_hits[-1].Start[0]/10.-offset[0], ar_muon_hits[-1].Start[1]/10.-offset[1], ar_muon_hits[-1].Start[2]/10.-offset[2])
            #print "Start z: ",hArStart.z()
            #print "End z: ",hArEnd.z()
            ar_dh = (hArEnd-hArStart).Mag()
            ar_dz = hArEnd.z()-hArStart.z()
                
            #print "ar_dh ", ar_dh
            #print "ar_dz ", ar_dz
            #if ar_dz == 0: continue
            ar_trk_length_gcm2 = 1.4 * ar_dh
            #print "Event # ", ient
            #print "AR Start Z ", hArStart.z()
            #print "AR End Z ", hArEnd.z()
            #print "AR TRACK LENGTH", ar_trk_length_gcm2

            

            #-------------------------------------------------------
            # look for muon hits in the scintillator
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "rmmsvol":
                    hits += key.second

            muon_hits = []
            for idx, hit in enumerate(hits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                if traj.ParentId == -1 and abs(traj.PDGCode) == 13:
                    muon_hits.append(hit)


            if len(muon_hits) < 2: continue

            hMuonStart = muon_hits[0].Start
            t_muonBirth[0] = hMuonStart[0]/10.-offset[0]
            t_muonBirth[1] = hMuonStart[1]/10.-offset[1]
            t_muonBirth[2] = hMuonStart[2]/10.-offset[2]

            trk_length_gcm2 = 0.
            de = 0.

            hPrev = None
            hFinal = None
            total_planes = 0

            for idx, hit in enumerate(muon_hits):
                de += hit.EnergyDeposit

                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )

                xpt.push_back(hStart.x())
                zpt.push_back(hStart.z())
                #print "hStart.z(): ", hStart.z()
                #print "hStart.x(): ", hStart.x()
                # this isn't the first hit, so we can start to build a track
                if hPrev is not None:
                    # thin layer is 3cm air + 1cm scint + 1.5cm steel = 5.5cm pitch
                    # thick layer is 3cm air + 1cm scint + 4cm steel = 8cm pitch
                    # if "gap" isn't a multiple of that, then something has gone off the rails and we want to stop
                    gap = hStart.z() - hPrev.z()
                    # correct for cosine of the incident angle
                    dh = (hStart-hPrev).Mag()
                    dz = hStart.z()-hPrev.z()
                    
                    # sometimes there are multiple hits within a scintillator plane, and we want to skip them
                    if gap < 4.:
                        continue
                    elif hStart.z() < 949.: # this is the "thin layer" region, so gap should be 5.5
                        nplanes = int(gap/5.5)
                        # if nplanes isn't an integer, then we've reached the end of the track
                        # and there is some other hit that is super far away typically that is messing up the length
                        if abs(nplanes - gap/5.5) > 1.E-6: 
                            break
                        if nplanes > nplanes_cut: # find an intolerable gap?
                            break
                        total_planes += nplanes
                        # we must have gone through 1cm scint + 1.5cm steel + 3cm air, divide by cos(theta)
                        trk_length_gcm2 += (1.05 + 1.5*7.85) * nplanes * dh/dz
                    elif (hStart.z() - 949.3) < 0.1: # the transition module, the gap is 9.5
                        if abs(gap-9.5) < 1.E-6:
                            trk_length_gcm2 += (1.05 + 4.*7.85) * dh/dz
                            total_planes += 1
                    else: #"thick layer region, gap should be 8cm
                        nplanes = int(gap/8.)
                        if abs(nplanes - gap/8.) > 1.E-6: # the track is gone to shit
                            break
                        if nplanes > nplanes_cut:
                            break
                        total_planes += nplanes

                        trk_length_gcm2 += (1.05 + 4.*7.85) * nplanes * dh/dz
                hPrev = hStart # update it
                hFinal = hStop

            t_muonDeath[0] = hFinal.x()
            t_muonDeath[1] = hFinal.y()
            t_muonDeath[2] = hFinal.z()

            t_muScintLen[0] = trk_length_gcm2 + extra_trk_length
            t_muLArLen[0] = ar_trk_length_gcm2 
            t_muScintEnergy[0] = de
            #print "Length ", t_muScintLen[0]
            #print "rmms KE ", t_rmmsKE[0]
            #print "exit KE ", t_muonExitKE[0]
            
            #ratio = 0. if not trk_length_gcm2 else t_rmmsKE[0]/trk_length_gcm2
            #print "%d,  %s Corrected Track Length %1.1f, out of Total muon KE %1.1f, (or ssri KE: %1.1f) planes %d, rmms length %1.1f, ar_trk_length %1.3f" % (ient, endVolName, 2.*t_muScintLen[0], t_lepKE[0], t_muonExitKE[0], total_planes, trk_length_gcm2, ar_trk_length_gcm2)

            tout.Fill()
        ient += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)

    (args, dummy) = parser.parse_args()

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
    t_lepRMMS = array('f',3*[0.0])
    tout.Branch('lepRMMS',t_lepRMMS,'leRMMS[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_rmmsKE = array('f',[0])
    tout.Branch('rmmsKE',t_rmmsKE,'rmmsKE/F')
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
    tout.Branch('muLArLen',t_muScintLen,'muLArLen/F')
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
    xpt = ROOT.std.vector('float')()
    zpt = ROOT.std.vector('float')()
    tout.Branch('xpt', xpt)
    tout.Branch('zpt', zpt)
    
    

    loaded = False
    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/%s.%d.edepsim.root" % (args.topdir, neutrino, run)
        print fname

        #fname = fname.replace("/pnfs","root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr")
        # see if it is an OK file
        #if not os.access( fname, os.R_OK ):
        #    print "Can't access file: %s" % fname
        #    continue
        tf = ROOT.TFile.Open( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            print "File is crap: %s" % fname
            continue

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )
        dspt.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    loop( events, dspt, tgeo, tout )

    fout.cd()
    tout.Write()
