import ROOT
from math import cos, atan2

geometries = ["rmms"]#["2Coils", "4Coils", "5Coils", "SPY"]
#topdir = "/dune/data/users/gsdavies/"
topdir="/dune/data/users/gsdavies/ssri/"
recs = ["contained", "rmms", "none"]

def loop( tree, hNum, hNumLen, hDenom, hNumX, hDenomX, hNumQ2, hDenomQ2, hVXY, hVYZ, hMuonBirthXZ, hMuonBirthYZ, hMuonDeathXZ, hMuonDeathYZ, hAngle, hED, hED2, hEDTrk, hScint, hMomAC, hMomSSRI, hMomCON, hMomRes, hMomRes2 ):

    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("vtx", 1)
    tree.SetBranchStatus("Ev", 1)
    tree.SetBranchStatus("lepPdg", 1)
    tree.SetBranchStatus("lepKE", 1)
    tree.SetBranchStatus("p3lep", 1)
    tree.SetBranchStatus("lepDeath", 1)
    tree.SetBranchStatus("muonDeath", 1)
    tree.SetBranchStatus("muonBirth", 1)
    tree.SetBranchStatus("muScintLen", 1)
    tree.SetBranchStatus("muScintEnergy", 1)
    tree.SetBranchStatus("muonReco", 1)
    tree.SetBranchStatus("muScintLenOld", 1)
    tree.SetBranchStatus("muonExitPt", 1)
    tree.SetBranchStatus("muonExitMom", 1)
    tree.SetBranchStatus("lepRMMS", 1)
    tree.SetBranchStatus("muonExitKE", 1)
    tree.SetBranchStatus("rmmsKE", 1)

    beam_angle = ROOT.TVector3( 0., 0., 1. )
    beam_angle.RotateX( 0.101 )

    N = tree.GetEntries()
    n = 0
    nevents = 0
    nmatch = 0
    nnu = 0
    nantinu = 0
    for entry in tree:
        if n % 1000 == 0:
            print "  event %d of %d..." % (n,N)
        n += 1

        # True numu CC
        if abs(entry.lepPdg) != 13: continue

        # vertex cut
        if abs(entry.vtx[0]) > 300. or abs(entry.vtx[1]) > 100. or entry.vtx[2] < 50. or entry.vtx[2] > 350.: continue

        #if (entry.lepPdg) == 13: continue
        mu = ROOT.TVector3( entry.p3lep[0], entry.p3lep[1], entry.p3lep[2] )
        muExit = ROOT.TVector3( entry.muonExitMom[0], entry.muonExitMom[1], entry.muonExitMom[2] )
        theta = mu.Angle(beam_angle)
        thetadeg = (180./3.1416)*theta


        angleA = 0.
        angleB = 0.
        diff   = 0.
        #vert = ROOT.TVector3 ( entry.vtx[0], entry.vtx[1], entry.vtx[2])
        #exitV = ROOT.TVector3 ( entry.muonExitPt[0], entry.muonExitPt[1], entry.muonExitPt[2] )
        #angleTanA = atan2(exitV.X(),exitV.Z()) - atan2(vert.X(),vert.Z())
        #print angleTanA
        deflection = entry.muonBirth[0] - entry.muonDeath[0]
        
        Elep = entry.lepKE/1000. + 0.105658
        plep = mu.Mag()/1000.
        plepexit = muExit.Mag()/1000.
        Q2 = 2*entry.Ev*(Elep - plep*cos(theta)) - 0.105658**2

        
        #print "plep: ", plep
        #print "plepexit: ", plepexit
        #if thetadeg < 20:
        hDenom.Fill( entry.lepKE/1000., thetadeg )
        hDenomQ2.Fill( entry.Ev, Q2 )

        #hVXY.Fill( entry.vtx[0], entry.vtx[1] )
        #hVYZ.Fill( entry.vtx[2], entry.vtx[1] )
       
        hMomAC.Fill(plep)


        rec = False
        if entry.muonReco == 1:
            hNum[0].Fill( entry.lepKE/1000., thetadeg )
            hNumLen[0].Fill( entry.lepKE/1000., entry.muScintLen )
            hNumQ2[0].Fill( entry.Ev, Q2 )
            hMomCON.Fill( plep )
            rec = True
        elif entry.muonReco == 0:
            hNum[2].Fill( entry.lepKE/1000., thetadeg )
            hNumQ2[2].Fill( entry.Ev, Q2 )
            rec = False
        elif entry.muonReco == 2:
            hNum[1].Fill( entry.lepKE/1000., thetadeg )
            #hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
            hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
            hNumQ2[1].Fill( entry.Ev, Q2 )
            rec = True
            # x -340
        match = False
        if entry.muonBirth[2] < 733. and entry.muonDeath[2] < 1300. and entry.muonDeath[1] < 25. and entry.muonDeath[1] > -260. and entry.muonDeath[0] < 300. and entry.muonDeath[0] > -300. and entry.muonExitKE > 0. and entry.lepPdg == 13:
            #if entry.muonReco == 2 and entry.lepPdg == -13:
            if (entry.muonDeath[0] < 165. and entry.muonDeath[0] > 10. or entry.muonDeath[0] > -165. and entry.muonDeath[0] < -10.):# and (entry.muonBirth[0] < 170. and entry.muonBirth[0] > 10. or entry.muonBirth > -170. and entry.muonBirth < -10.):

                #print "Deflection: ", deflection
                #print "PDG: ", entry.lepPdg

                if deflection < 0:
                    if entry.lepPdg == 13:
                        match = True
                    elif entry.lepPdg == -13:
                        match = False
                elif deflection > 0:
                    if entry.lepPdg == 13:
                        match = False
                    elif entry.lepPdg == -13:
                        match = True
                        
                nevents += 1
                if match == True: nmatch += 1
                if entry.lepPdg == 13: nnu += 1
                if entry.lepPdg == -13: nantinu += 1
            
                hVXY.Fill( entry.lepDeath[0], entry.lepDeath[1] )
                hVYZ.Fill( entry.lepDeath[2], entry.lepDeath[1] )
                #if entry.rmmsKE/1000. > 2. and entry.muScintLen < 500.:
                hMuonDeathYZ.Fill( entry.muonDeath[2], entry.muonDeath[1] )
                hMuonDeathXZ.Fill( entry.muonDeath[2], entry.muonDeath[0] )
                hMuonBirthYZ.Fill( entry.muonBirth[2], entry.muonBirth[1] )
                hMuonBirthXZ.Fill( entry.muonBirth[2], entry.muonBirth[0] )
            
                hED.Fill( entry.muonExitKE/1000., entry.muScintLen )
                #reconstructed mom vs. "true" 
                hED2.Fill( plepexit, entry.muScintLen )
                hMomRes2.Fill(plepexit, ((entry.muScintLen+5.42)/4963.)-plepexit)
                #if plepexit > 1. and plepexit < 1.5:
                hMomRes.Fill( entry.muonExitKE/1000., (2*entry.muScintLen-entry.muonExitKE)/entry.muonExitKE )
                dedx = 0.
                if entry.muScintLen == 0.:
                    dedx = 0.
                else: dedx = entry.muScintEnergy/entry.muScintLen
                hEDTrk.Fill(dedx, plepexit )
                hScint.Fill(entry.muScintEnergy)
                hAngle.Fill( entry.muScintLen, diff )
                hMomSSRI.Fill(plepexit)
                rec = True
        print "Nevents: ", nevents
        print "N match: ", nmatch
        print "N anti-nu: ", nantinu
        print "N Nu: ", nnu

        if thetadeg < 30.:
            hDenomX.Fill( entry.lepKE/1000., entry.vtx[0] )
            if rec:
                hNumX.Fill( entry.lepKE/1000., entry.vtx[0] )

if __name__ == "__main__":

    hDenom = [None for g in geometries]

    hNum = [[None for r in recs] for g in geometries]
    hNumLen = [[None for r in recs] for g in geometries]
    hDenomQ2 = [None for g in geometries]
    hNumQ2 = [[None for r in recs] for g in geometries]
    hDenomX = [None for g in geometries]
    hNumX = [None for g in geometries]
    hVXY = [None for g in geometries]
    hVYZ = [None for g in geometries]
    #hVXZ = [None for g in geometries]
    hMuonBirthXZ = [None for g in geometries]
    hMuonBirthYZ = [None for g in geometries]
    hMuonDeathXZ = [None for g in geometries]
    hMuonDeathYZ = [None for g in geometries]
    hAngle = [None for g in geometries]
    hED = [None for g in geometries]
    hED2 = [None for g in geometries]
    hEDTrk = [None for g in geometries]
    hScint = [None for g in geometries]
    hMomAC = [None for g in geometries]
    hMomSSRI = [None for g in geometries]
    hMomCON = [None for g in geometries]
    hMomRes = [None for g in geometries]
    hMomRes2 = [None for g in geometries]
    for g,geom in enumerate(geometries):
        tree = ROOT.TChain("tree","tree")
        #tree.Add( "%s/MPD_%s_LAr/*.root" % (topdir, geom) )
        #tree.Add("%s/neutrino.combined.dump.root" % topdir )
        #tree.Add("%s/combined_fifty.root" % topdir )
        #tree.Add("%s/dump_uniform_10.root" % topdir )
        tree.Add("%s/dump_cmfix_allnu-v1.root" % topdir )
        #tree.Add("%s/dump_cmfix_allnu.root" % topdir )

        hAngle[g] = ROOT.TH2D( "angle_%s" % geom, ";Neutrino energy (GeV);Angle (RMMS) [degrees]", 200, 0., 2000., 45, -90., 90. )
        hED[g] = ROOT.TH2D( "MOM-EDEP_%s" % geom, ";Muon Kinetic Energy (GeV) [SSRI];Muon Track Length (g/cm^2)", 35, 0., 7., 40, 0., 4000. )
        hED2[g] = ROOT.TH2D( "MOMorig-EDEP_%s" % geom, ";Muon Momentum (GeV [SSRI] ;Reconstructed Muon Momentum (GeV) [SSRI]", 60, 0., 6., 100, 0., 2000. )
        hEDTrk[g] = ROOT.TH2D( "edep-trk_%s" % geom, ";Muon Momentum (GeV) [SSRI];#mu dE/dx (MeV/g/cm^{2})", 600, 0.05, 2., 500, 0.0, 5. )
        hScint[g] = ROOT.TH1D( "scint_%s" % geom, ";Energy Deposited", 100, 0., 500. )
        hMomAC[g] = ROOT.TH1D( "momAC_%s" % geom, ";Muon Momentum (All)", 200, 0., 20. )
        hMomSSRI[g] = ROOT.TH1D( "momSSRI_%s" % geom, ";Muon Momentum (SSRI)", 200, 0., 20. )
        hMomCON[g] = ROOT.TH1D( "momCON_%s" % geom, ";Muon Momentum (ArgonCube)", 200, 0., 20. )
        hMomRes2[g] = ROOT.TH2D( "momRES2_%s" % geom, ";True Muon Momentum (GeV) [SSRI];(Reco - True) Muom Mom. (GeV) [SSRI]",  20, 0., 10., 5000, -10., 10.)
        #hMomRes[g] = ROOT.TH2D( "momRES_%s" % geom, ";True Momentum (SSRI);Momentum Resolution (SSRI)", 40, 0., 20., 100, 0., 1. )
        hMomRes[g] = ROOT.TH2D( "momRES_%s" % geom, ";#mu KE (GeV) [SSRI];Fractional Residual", 20, 0., 5., 100, -1.,1. )

        hVXY[g] = ROOT.TH2D( "denomvxy_%s" % geom, ";Vertex X;vertex Y", 40, -350., 350., 40, -400., 400. )
        hVYZ[g] = ROOT.TH2D( "denomvyz_%s" % geom, ";Vertex Z;vertex Y", 100, 0., 2000., 40, -400., 400. )
        hMuonBirthYZ[g] = ROOT.TH2D( "birthyz_%s" % geom, ";#mu Entry Z (cm); #mu Entry Y (cm)", 100, 500., 1500., 30, -400., 200. )
        hMuonBirthXZ[g] = ROOT.TH2D( "birthxz_%s" % geom, ";#mu Entry Z (cm);#mu Entry X (cm)", 100, 500., 1500., 40, -400., 400. )
        hMuonDeathYZ[g] = ROOT.TH2D( "deathyz_%s" % geom, ";#mu Endpoint Z (cm); #mu Endpoint Y (cm)", 100, 500., 1500., 30, -400., 200. )
        hMuonDeathXZ[g] = ROOT.TH2D( "deathxz_%s" % geom, ";#mu Endpoint Z (cm);#mu Endpoint X (cm)", 100, 500., 1500., 1600, -400., 400. )
        hDenom[g] = ROOT.TH2D( "denom_%s" % geom, ";Muon kinetic energy (GeV);Muon angle (degrees)", 30, 0., 6., 18, 0., 180. )
        hDenomQ2[g] = ROOT.TH2D( "denomQ2_%s" % geom, ";Neutrino energy (GeV);Q^{2} (GeV^{2})", 25, 0., 5., 25, 0., 5. )

        hDenomX[g] = ROOT.TH2D( "denomX_%s" % geom, ";Muon kinetic energy (GeV);Vertex X (cm)", 30, 0., 6., 60, -300., 300. )
        hNumX[g] = ROOT.TH2D( "numX_%s" % geom, ";Muon kinetic energy (GeV);Vertex X (cm)", 30, 0., 6., 60, -300., 300. )
        for r,rec in enumerate(recs):
            hNum[g][r] = ROOT.TH2D( "num_%s_%s" % (geom, rec), ";Muon kinetic energy (GeV);Muon angle (degrees)", 30, 0., 6., 18, 0., 180. )
            hNumLen[g][r] = ROOT.TH2D( "numlen_%s_%s" % (geom, rec), ";Muon kinetic energy (GeV);Muon track length", 30, 0., 6., 20, 0., 2000. )
            hNumQ2[g][r] = ROOT.TH2D( "numQ2_%s_%s" % (geom, rec), ";Neutrino energy (GeV);Q^{2} (GeV^{2})", 25, 0., 5., 25, 0., 5. )
        print "Looping for %s" % geom
        loop( tree, hNum[g], hNumLen[g], hDenom[g], hNumX[g], hDenomX[g], hNumQ2[g], hDenomQ2[g], hVXY[g], hVYZ[g], hMuonBirthXZ[g], hMuonBirthYZ[g], hMuonDeathXZ[g], hMuonDeathYZ[g], hAngle[g], hED[g], hED2[g], hEDTrk[g], hScint[g], hMomAC[g], hMomSSRI[g], hMomCON[g], hMomRes[g], hMomRes2[g] )

        #out = ROOT.TFile( "%s/filled_fifty.root" % topdir, "RECREATE" )
    out = ROOT.TFile( "%s/dumpall_cmfixhistos.root" % topdir, "RECREATE" )
    for g,geom in enumerate(geometries):
        hDenom[g].Write()
        hDenomQ2[g].Write()
        hNumX[g].Write()
        hDenomX[g].Write()
        hVXY[g].Write()
        hVYZ[g].Write()
        hAngle[g].Write()
        hED[g].Write()
        hED2[g].Write()
        hEDTrk[g].Write()
        hScint[g].Write()
        hMomAC[g].Write()
        hMomSSRI[g].Write()
        hMomCON[g].Write()
        hMomRes[g].Write()
        hMomRes2[g].Write()
        hMuonBirthXZ[g].Write()
        hMuonBirthYZ[g].Write()
        hMuonDeathXZ[g].Write()
        hMuonDeathYZ[g].Write()

        for r,rec in enumerate(recs):
            hNum[g][r].Write()
            hNumLen[g][r].Write()
            hNumQ2[g][r].Write()


