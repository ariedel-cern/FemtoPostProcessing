/*
 * File              : postProcessing.C
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 24.08.2022
 * Last Modified Date: 10.11.2022
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

#include <RtypesCore.h>
#include <TDataType.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

Float_t ConvertBin(Char_t Input);

Int_t postProcessing(const char *ConfigFile, const char *DataFile,
                     const char *OutputFile) {

  // load config file
  std::fstream Jfile(ConfigFile);
  nlohmann::json Jconfig = nlohmann::json::parse(Jfile);

  // sane defaults for histograms
  Float_t ptRangeMin = Jconfig["ptRangeMin"].get<Float_t>();
  Float_t ptRangeMax = Jconfig["ptRangeMax"].get<Float_t>();
  Float_t etaRangeMin = -1., etaRangeMax = 1., phiRangeMin = 0.,
          phiRangeMax = TMath::TwoPi(), dcazRangeMin = -0.3, dcazRangeMax = 0.3,
          dcaxyRangeMin = -0.3, dcaxyRangeMax = 0.3, daughdcaRangeMin = -0.3,
          daughdcaRangeMax = 0.3, nsigmaTPCRangeMin = -8.,
          transradiusRangeMin = 0, transradiusRangeMax = 150,
          nsigmaTPCRangeMax = 8., nsigmaTOFRangeMin = -8.,
          nsigmaTOFRangeMax = 8., tpcsignalRangeMin = 0,
          tpcsignalRangeMax = 500, invMassLambda0 = 0., invMassLambda1 = 2.;

  Int_t ptBins = Jconfig["ptBins"].get<Float_t>();
  Int_t etaBins = 1000, phiBins = 1000, dcazBins = 300, dcaxyBins = 300,
        daughdcaBins = 300, transradiusBins = 1000, nsigmaTPCBins = 100,
        nsigmaTOFBins = 100, tpcsignalBins = 500, invMassBins = 100;

  TH1F *ptDeuteron =
      new TH1F("ptDeuteron", "ptDeuteron", ptBins, ptRangeMin, ptRangeMax);
  TH1F *phiDeuteron =
      new TH1F("phiDeuteron", "phiDeuteron", phiBins, phiRangeMin, phiRangeMax);
  TH1F *etaDeuteron =
      new TH1F("etaDeuteron", "etaDeuteron", etaBins, etaRangeMin, etaRangeMax);
  TH1F *dcazDeuteron = new TH1F("dcazDeuteron", "dcazDeuteron", dcazBins,
                                dcazRangeMin, dcazRangeMax);
  TH1F *dcaxyDeuteron = new TH1F("dcaxyDeuteron", "dcaxyDeuteron", dcaxyBins,
                                 dcaxyRangeMin, dcaxyRangeMax);
  TH2F *dcaz_pt_Deuteron =
      new TH2F("dcaz_pt_Deuteron", "dcaz_pt_Deuteron", ptBins, ptRangeMin,
               ptRangeMax, dcazBins, dcazRangeMin, dcazRangeMax);
  TH2F *dcaxy_pt_Deuteron =
      new TH2F("dcaxy_pt_Deuteron", "dcaxy_pt_Deuteron", ptBins, ptRangeMin,
               ptRangeMax, dcaxyBins, dcaxyRangeMin, dcaxyRangeMax);
  TH1F *tpcDeuteron =
      new TH1F("nsigmatpcDeuteron", "nsigmatpcDeuteron", nsigmaTPCBins,
               nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH1F *tofDeuteron =
      new TH1F("nsigmatofDeuteron", "nsigmatofDeuteron", nsigmaTOFBins,
               nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH1F *tpcsignalDeuteron =
      new TH1F("tpcsignalDeuteron", "tpcsignalDeuteron", tpcsignalBins,
               tpcsignalRangeMin, tpcsignalRangeMax);
  TH2F *tpc_pt_Deuteron = new TH2F(
      "nsigmatpc_pt_Deuteron", "nsigmatpc_pt_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, nsigmaTPCBins, nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH2F *tof_pt_Deuteron = new TH2F(
      "nsigmatof_pt_Deuteron", "nsigmatof_pt_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, nsigmaTOFBins, nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH2F *tpcsignal_pt_Deuteron = new TH2F(
      "tpcsignal_pt_Deuteron", "tpcsignal_pt_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, tpcsignalBins, tpcsignalRangeMin, tpcsignalRangeMax);
  TH2F *tpc_p_Deuteron = new TH2F(
      "nsigmatpc_p_Deuteron", "nsigmatpc_p_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, nsigmaTPCBins, nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH2F *tof_p_Deuteron = new TH2F(
      "nsigmatof_p_Deuteron", "nsigmatof_p_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, nsigmaTOFBins, nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH2F *tpcsignal_p_Deuteron = new TH2F(
      "tpcsignal_p_Deuteron", "tpcsignal_p_Deuteron", ptBins, ptRangeMin,
      ptRangeMax, tpcsignalBins, tpcsignalRangeMin, tpcsignalRangeMax);

  TH1F *ptProton =
      new TH1F("ptProton", "ptProton", ptBins, ptRangeMin, ptRangeMax);
  TH1F *phiProton =
      new TH1F("phiProton", "phiProton", phiBins, phiRangeMin, phiRangeMax);
  TH1F *etaProton =
      new TH1F("etaProton", "etaProton", etaBins, etaRangeMin, etaRangeMax);
  TH1F *dcazProton = new TH1F("dcazProton", "dcazProton", dcazBins,
                              dcazRangeMin, dcazRangeMax);
  TH1F *dcaxyProton = new TH1F("dcaxyProton", "dcaxyProton", dcaxyBins,
                               dcaxyRangeMin, dcaxyRangeMax);
  TH2F *dcaz_pt_Proton =
      new TH2F("dcaz_pt_Proton", "dcaz_pt_Proton", ptBins, ptRangeMin,
               ptRangeMax, dcazBins, dcazRangeMin, dcazRangeMax);
  TH2F *dcaxy_pt_Proton =
      new TH2F("dcaxy_pt_Proton", "dcaxy_pt_Proton", ptBins, ptRangeMin,
               ptRangeMax, dcaxyBins, dcaxyRangeMin, dcaxyRangeMax);
  TH1F *tpcProton =
      new TH1F("nsigmatpcProton", "nsigmatpcProton", nsigmaTPCBins,
               nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH1F *tofProton =
      new TH1F("nsigmatofProton", "nsigmatofProton", nsigmaTOFBins,
               nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH1F *tpcsignalProton =
      new TH1F("tpcsignalProton", "tpcsignalProton", tpcsignalBins,
               tpcsignalRangeMin, tpcsignalRangeMax);
  TH2F *tpc_pt_Proton =
      new TH2F("nsigmatpc_pt_Proton", "nsigmatpc_pt_Proton", ptBins, ptRangeMin,
               ptRangeMax, nsigmaTPCBins, nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH2F *tof_pt_Proton =
      new TH2F("nsigmatof_pt_Proton", "nsigmatof_pt_Proton", ptBins, ptRangeMin,
               ptRangeMax, nsigmaTOFBins, nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH2F *tpcsignal_pt_Proton =
      new TH2F("tpcsignal_pt_Proton", "tpcsignal_pt_Proton", ptBins, ptRangeMin,
               ptRangeMax, tpcsignalBins, tpcsignalRangeMin, tpcsignalRangeMax);
  TH2F *tpc_p_Proton =
      new TH2F("nsigmatpc_p_Proton", "nsigmatpc_p_Proton", ptBins, ptRangeMin,
               ptRangeMax, nsigmaTPCBins, nsigmaTPCRangeMin, nsigmaTPCRangeMax);
  TH2F *tof_p_Proton =
      new TH2F("nsigmatof_p_Proton", "nsigmatof_p_Proton", ptBins, ptRangeMin,
               ptRangeMax, nsigmaTOFBins, nsigmaTOFRangeMin, nsigmaTOFRangeMax);
  TH2F *tpcsignal_p_Proton =
      new TH2F("tpcsignal_p_Proton", "tpcsignal_p_Proton", ptBins, ptRangeMin,
               ptRangeMax, tpcsignalBins, tpcsignalRangeMin, tpcsignalRangeMax);

  TH1F *ptLambda =
      new TH1F("ptLambda", "ptLambda", ptBins, ptRangeMin, ptRangeMax);
  TH1F *phiLambda =
      new TH1F("phiLambda", "phiLambda", phiBins, phiRangeMin, phiRangeMax);
  TH1F *etaLambda =
      new TH1F("etaLambda", "etaLambda", etaBins, etaRangeMin, etaRangeMax);
  TH1F *daughDCALambda =
      new TH1F("daughDCALambda", "daughDCALambda", daughdcaBins,
               daughdcaRangeMin, daughdcaRangeMax);
  TH1F *transradiusLambda =
      new TH1F("transradiusLambda", "transradiusLambda", transradiusBins,
               transradiusRangeMin, transradiusRangeMax);
  TH1F *HistInvMassLambda =
      new TH1F("invMassLambda", "invMassLambda", invMassBins, invMassLambda0,
               invMassLambda1);

  TH1F *ptPosDaugh =
      new TH1F("ptPosDaugh", "ptPosDaugh", ptBins, ptRangeMin, ptRangeMax);
  TH1F *phiPosDaugh =
      new TH1F("phiPosDaugh", "phiPosDaugh", phiBins, phiRangeMin, phiRangeMax);
  TH1F *etaPosDaugh =
      new TH1F("etaPosDaugh", "etaPosDaugh", etaBins, etaRangeMin, etaRangeMax);
  TH1F *nsigmaTPCPosDaugh =
      new TH1F("nsigmaTPCPosDaugh", "nsigmaTPCPosDaugh", nsigmaTPCBins,
               nsigmaTPCRangeMin, nsigmaTPCRangeMax);

  TH1F *ptNegDaugh =
      new TH1F("ptNegDaugh", "ptNegDaugh", ptBins, ptRangeMin, ptRangeMax);
  TH1F *phiNegDaugh =
      new TH1F("phiNegDaugh", "phiNegDaugh", phiBins, phiRangeMin, phiRangeMax);
  TH1F *etaNegDaugh =
      new TH1F("etaNegDaugh", "etaNegDaugh", etaBins, etaRangeMin, etaRangeMax);
  TH1F *nsigmaTPCNegDaugh =
      new TH1F("nsigmaTPCNegDaugh", "nsigmaTPCNegDaugh", nsigmaTPCBins,
               nsigmaTPCRangeMin, nsigmaTPCRangeMax);

  TH1F *HistDaughterPt =
      new TH1F("ptDaughter", "ptDaughter", ptBins, ptRangeMin, ptRangeMax);

  Bool_t UsePid = Jconfig["UsePid"].get<Bool_t>();
  Float_t posZMax = Jconfig["posZMax"].get<Float_t>();
  Float_t etaMax = Jconfig["etaMax"].get<Float_t>();
  Float_t dcazMax = Jconfig["dcazMax"].get<Float_t>();
  Float_t dcaxyMax = Jconfig["dcaxyMax"].get<Float_t>();
  Float_t TPCclustersMin = Jconfig["TPCclustersMin"].get<Float_t>();
  Float_t TPCcrossedrowsMin = Jconfig["TPCcrossedrowsMin"].get<Float_t>();
  Float_t TPCcrossedrowsOverclustersMin =
      Jconfig["TPCcrossedrowsOverclustersMin"].get<Float_t>();
  Float_t ITSclustersMin = Jconfig["ITSclustersMin"].get<Float_t>();
  Float_t ITSclustersIBMin = Jconfig["ITSclustersIBMin"].get<Float_t>();
  Float_t nsigmaTPCDeuteronMax = Jconfig["nsigmaTPCDeuteronMax"].get<Float_t>();
  Float_t nsigmaTPCRejection = Jconfig["nsigmaTPCRejection"].get<Float_t>();
  Float_t nsigmaTPCProtonMax = Jconfig["nsigmaTPCProtonMax"].get<Float_t>();
  Float_t nsigmaTPCTOFProtonMax =
      Jconfig["nsigmaTPCTOFProtonMax"].get<Float_t>();
  Float_t pPIDThresholdProton = Jconfig["pPIDThresholdProton"].get<Float_t>();
  Float_t ptDeuteronMin = Jconfig["ptDeuteronMin"].get<Float_t>();
  Float_t ptDeuteronMax = Jconfig["ptDeuteronMax"].get<Float_t>();
  Float_t ptProtonMin = Jconfig["ptProtonMin"].get<Float_t>();
  Float_t ptProtonMax = Jconfig["ptProtonMax"].get<Float_t>();
  Float_t ptLambdaMin = Jconfig["ptLambdaMin"].get<Float_t>();
  Float_t ptLambdaMax = Jconfig["ptLambdaMax"].get<Float_t>();
  Float_t daughDCAMax = Jconfig["daughDCAMax"].get<Float_t>();
  Float_t TransRadiusMin = Jconfig["TransRadiusMin"].get<Float_t>();
  Float_t TransRadiusMax = Jconfig["TransRadiusMax"].get<Float_t>();

  // load data file
  TFile *file = new TFile(DataFile, "READ");
  TDirectoryFile *TDirFile;
  TTree *TreeDebugParts, *TreePosDaughDebug, *TreeNegDaughDebug, *TreeParts,
      *TreePosDaugh, *TreeNegDaugh, *TreeCols;

  // tree variables
  Float_t pt, p, phi, eta, dcaz, dcaxy, daughDCA, DecayVtxX, DecayVtxY,
      DecayVtxZ, TransRadius, signalTPC, posZ, mult, mLambda, mAntiLambda,
      mKaon;
  Float_t PosDaugh_pt, PosDaugh_phi, PosDaugh_eta, PosDaugh_TPCnCls,
      PosDaugh_nsigmaTPC;
  Float_t NegDaugh_pt, NegDaugh_phi, NegDaugh_eta, NegDaugh_TPCnCls,
      NegDaugh_nsigmaTPC;

  uint8_t nTPCClusters, nTPCFindable, nTPCCrossedRows, nTPCClustersShared,
      nITSClusters, nITSClustersInnerBarrel, PartType;
  Char_t sign, nSigmaTPCDeuteron, nSigmaTPCProton, nSigmaTOFDeuteron,
      nSigmaTOFProton, nSigmaTPCElectron, nSigmaTPCPion, nSigmaTPCKaon;

  int32_t CollisionID;
  std::vector<uint32_t> passedCollisionIDs = {};

  bool keepProton = false;
  bool keepDeuteron = false;

  TH1F *histPosz = new TH1F("posz", "posz", 1000, -20, 20);
  TH1F *histMult = new TH1F("mul", "mul", 10000, 0, 10000);

  for (TObject *key : *(file->GetListOfKeys())) {
    TDirFile = dynamic_cast<TDirectoryFile *>(file->Get(key->GetName()));

    if (!TDirFile) {
      std::cout << "Did not get a valid TDirectoryFile. Skip..." << std::endl;
      continue;
    }

    std::cout << "Working on TDirFile " << TDirFile->GetName() << std::endl;

    TreeParts = dynamic_cast<TTree *>(TDirFile->Get("O2femtodreamparts"));
    TreePosDaugh = dynamic_cast<TTree *>(TDirFile->Get("O2femtodreamparts"));
    TreeNegDaugh = dynamic_cast<TTree *>(TDirFile->Get("O2femtodreamparts"));
    TreeCols = dynamic_cast<TTree *>(TDirFile->Get("O2femtodreamcols"));
    TreeDebugParts = dynamic_cast<TTree *>(TDirFile->Get("O2femtodebugparts"));

    if (!TreeParts || !TreeDebugParts || !TreeCols) {
      std::cout << "One of the trees was not found. Skip..." << std::endl;
      continue;
    }

    TreeCols->SetBranchAddress("fPosZ", &posZ);
    TreeCols->SetBranchAddress("fMultV0M", &mult);

    TreeParts->SetBranchAddress("fPt", &pt);
    TreeParts->SetBranchAddress("fPhi", &phi);
    TreeParts->SetBranchAddress("fEta", &eta);
    TreeParts->SetBranchAddress("fIndexFemtoDreamCollisions", &CollisionID);
    TreeParts->SetBranchAddress("fPartType", &PartType);
    TreeParts->SetBranchAddress("fMLambda", &mLambda);
    TreeParts->SetBranchAddress("fMAntiLambda", &mAntiLambda);

    TreePosDaugh->SetBranchAddress("fPt", &PosDaugh_pt);
    TreePosDaugh->SetBranchAddress("fPhi", &PosDaugh_phi);
    TreePosDaugh->SetBranchAddress("fEta", &PosDaugh_eta);

    TreeNegDaugh->SetBranchAddress("fPt", &NegDaugh_pt);
    TreeNegDaugh->SetBranchAddress("fPhi", &NegDaugh_phi);
    TreeNegDaugh->SetBranchAddress("fEta", &NegDaugh_eta);

    TreeDebugParts->SetBranchAddress("fSign", &sign);
    TreeDebugParts->SetBranchAddress("fDcaZ", &dcaz);
    TreeDebugParts->SetBranchAddress("fDcaXY", &dcaxy);
    TreeDebugParts->SetBranchAddress("fDaughDCA", &daughDCA);
    TreeDebugParts->SetBranchAddress("fDecayVtxX", &DecayVtxX);
    TreeDebugParts->SetBranchAddress("fDecayVtxY", &DecayVtxY);
    TreeDebugParts->SetBranchAddress("fDecayVtxZ", &DecayVtxZ);
    TreeDebugParts->SetBranchAddress("fTransRadius", &TransRadius);
    TreeDebugParts->SetBranchAddress("fMKaon", &mKaon);
    TreeDebugParts->SetBranchAddress("fITSNCls", &nITSClusters);
    TreeDebugParts->SetBranchAddress("fITSNClsInnerBarrel",
                                     &nITSClustersInnerBarrel);
    TreeDebugParts->SetBranchAddress("fTPCNClsFound", &nTPCClusters);
    TreeDebugParts->SetBranchAddress("fTPCNClsFindable", &nTPCFindable);
    TreeDebugParts->SetBranchAddress("fTPCNClsShared", &nTPCClustersShared);
    TreeDebugParts->SetBranchAddress("fTPCNClsCrossedRows", &nTPCCrossedRows);
    TreeDebugParts->SetBranchAddress("fTPCNSigmaStoreEl", &nSigmaTPCElectron);
    TreeDebugParts->SetBranchAddress("fTPCNSigmaStorePi", &nSigmaTPCPion);
    TreeDebugParts->SetBranchAddress("fTPCNSigmaStoreKa", &nSigmaTPCKaon);
    TreeDebugParts->SetBranchAddress("fTPCNSigmaStoreDe", &nSigmaTPCDeuteron);
    TreeDebugParts->SetBranchAddress("fTOFNSigmaStoreDe", &nSigmaTOFDeuteron);
    TreeDebugParts->SetBranchAddress("fTPCNSigmaStorePr", &nSigmaTPCProton);
    TreeDebugParts->SetBranchAddress("fTOFNSigmaStorePr", &nSigmaTOFProton);
    TreeDebugParts->SetBranchAddress("fTPCSignal", &signalTPC);

    TreePosDaughDebug->SetBranchAddress("fTPCNClsFound", &PosDaugh_TPCnCls);
    TreePosDaughDebug->SetBranchAddress("fTPCNSigmaStorePr",
                                        &PosDaugh_nsigmaTPC);

    TreePosDaughDebug->SetBranchAddress("fTPCNClsFound", &NegDaugh_TPCnCls);
    TreePosDaughDebug->SetBranchAddress("fTPCNSigmaStorePi",
                                        &NegDaugh_nsigmaTPC);

    // loop over all particles
    for (int i = 0; i < TreeParts->GetEntries(); i++) {

      TreeParts->GetEntry(i);
      p = pt * std::cosh(eta);
      TreeDebugParts->GetEntry(i);

      // get corresponding collision
      TreeCols->GetEntry(CollisionID);

      if (std::abs(posZ) > posZMax) {
        continue;
      }

      // fill posz, but check if we filled the collision before
      if (std::find(passedCollisionIDs.begin(), passedCollisionIDs.end(),
                    CollisionID) == passedCollisionIDs.end()) {
        histPosz->Fill(posZ);
        histMult->Fill(mult);
        passedCollisionIDs.push_back(CollisionID);
      }

      // lambda
      if (PartType == static_cast<uint8_t>(1)) {

        // cuts on lambda
        if (sign < 0 || pt < ptLambdaMin || pt > ptLambdaMax ||
            daughDCA > daughDCAMax || TransRadius < TransRadiusMin ||
            TransRadius > TransRadiusMax || std::abs(eta) > etaMax) {
          continue;
        }

        // cut on lambda daughters
        TreePosDaugh->GetEntry(i + 1);
        TreePosDaughDebug->GetEntry(i + 1);
        TreeNegDaugh->GetEntry(i + 2);
        TreeNegDaughDebug->GetEntry(i + 2);

        if (std::abs(PosDaugh_nsigmaTPC) > 6 ||
            std::abs(NegDaugh_nsigmaTPC) > 6 ||
            PosDaugh_TPCnCls < DaughTPCnclsMin ||
            NegDaugh_TPCnCls < DaughTPCnclsMin) {
          continue;
        }

        ptLambda->Fill(pt);
        etaLambda->Fill(eta);
        phiLambda->Fill(phi);
        daughDCALambda->Fill(daughDCA);
        transradiusLambda->Fill(TransRadius);
        HistInvMassLambda->Fill(mLambda);

        continue;
      }

      // daughters already handled in lambda
      // if (PartType == static_cast<uint8_t>(2)) {
      //   HistDaughterPt->Fill(pt);
      //   continue;
      // }

      // if not a track, bail out immediately
      // process protons and deuterons after this point
      if (PartType != static_cast<uint8_t>(0)) {
        continue;
      }

      // general cuts for protons and deuterons
      if (static_cast<Int_t>(nTPCClusters) < TPCclustersMin ||
          static_cast<Int_t>(nTPCCrossedRows) < TPCcrossedrowsMin ||
          static_cast<Int_t>(nTPCClustersShared) > 0 ||
          std::abs(eta) > etaMax ||
          (1.f * static_cast<Int_t>(nTPCCrossedRows) /
           static_cast<Int_t>(nTPCFindable)) < TPCcrossedrowsOverclustersMin ||
          static_cast<Int_t>(sign) < 0) {
        continue;
      }

      keepProton = false;
      keepDeuteron = false;

      // deuteron
      if ((UsePid &&
           std::abs(ConvertBin(nSigmaTPCDeuteron)) < nsigmaTPCDeuteronMax) &&
          pt > ptDeuteronMin &&
          pt<ptDeuteronMax &&static_cast<Int_t>(nITSClusters)> ITSclustersMin &&
          static_cast<Int_t>(nITSClustersInnerBarrel) > ITSclustersIBMin) {

        if (nsigmaTPCRejection > 0) {
          if (std::abs(ConvertBin(nSigmaTPCElectron)) > nsigmaTPCRejection &&
              std::abs(ConvertBin(nSigmaTPCPion)) > nsigmaTPCRejection &&
              std::abs(ConvertBin(nSigmaTPCProton)) > nsigmaTPCRejection) {
            keepDeuteron = true;
          }
        } else {
          keepDeuteron = true;
        }
      }

      if (keepDeuteron) {
        ptDeuteron->Fill(pt);
        phiDeuteron->Fill(phi);
        etaDeuteron->Fill(eta);
        dcazDeuteron->Fill(dcaz);
        dcaxyDeuteron->Fill(dcaxy);
        dcaz_pt_Deuteron->Fill(pt, dcaz);
        dcaxy_pt_Deuteron->Fill(pt, dcaxy);
        tpcDeuteron->Fill(ConvertBin(nSigmaTPCDeuteron));
        tofDeuteron->Fill(ConvertBin(nSigmaTOFDeuteron));
        tpcsignalDeuteron->Fill(signalTPC);
        tpc_pt_Deuteron->Fill(pt, ConvertBin(nSigmaTPCDeuteron));
        tof_pt_Deuteron->Fill(pt, ConvertBin(nSigmaTOFDeuteron));
        tpcsignal_pt_Deuteron->Fill(pt, signalTPC);
        tpc_p_Deuteron->Fill(p, ConvertBin(nSigmaTPCDeuteron));
        tof_p_Deuteron->Fill(p, ConvertBin(nSigmaTOFDeuteron));
        tpcsignal_p_Deuteron->Fill(p, signalTPC);
      }

      // proton
      if (pt < ptProtonMax && pt > ptProtonMin) {
        if (p < pPIDThresholdProton) {
          if (UsePid &&
              std::abs(ConvertBin(nSigmaTPCProton)) < nsigmaTPCProtonMax) {
            keepProton = true;
          }
        } else {
          if (UsePid &&
              (TMath::Sqrt(TMath::Power(ConvertBin(nSigmaTPCProton), 2) +
                           TMath::Power(ConvertBin(nSigmaTOFProton), 2))) <
                  nsigmaTPCTOFProtonMax) {
            keepProton = true;
          }
        }
      }

      if (keepProton) {
        ptProton->Fill(pt);
        phiProton->Fill(phi);
        etaProton->Fill(eta);
        dcazProton->Fill(dcaz);
        dcaxyProton->Fill(dcaxy);
        dcaz_pt_Proton->Fill(pt, dcaz);
        dcaxy_pt_Proton->Fill(pt, dcaxy);
        tpcProton->Fill(ConvertBin(nSigmaTPCProton));
        tofProton->Fill(ConvertBin(nSigmaTOFProton));
        tpcsignalProton->Fill(signalTPC);
        tpc_pt_Proton->Fill(pt, ConvertBin(nSigmaTPCProton));
        tof_pt_Proton->Fill(pt, ConvertBin(nSigmaTOFProton));
        tpcsignal_pt_Proton->Fill(pt, signalTPC);
        tpc_p_Proton->Fill(p, ConvertBin(nSigmaTPCProton));
        tof_p_Proton->Fill(p, ConvertBin(nSigmaTOFProton));
        tpcsignal_p_Proton->Fill(p, signalTPC);
      }
    }
  }

  TFile *Output = new TFile(OutputFile, "RECREATE");

  TList *DeuteronList = new TList();

  DeuteronList->Add(ptDeuteron);
  DeuteronList->Add(phiDeuteron);
  DeuteronList->Add(etaDeuteron);
  DeuteronList->Add(dcazDeuteron);
  DeuteronList->Add(dcaxyDeuteron);
  DeuteronList->Add(dcaz_pt_Deuteron);
  DeuteronList->Add(dcaxy_pt_Deuteron);
  DeuteronList->Add(tpcDeuteron);
  DeuteronList->Add(tofDeuteron);
  DeuteronList->Add(tpcsignalDeuteron);
  DeuteronList->Add(tpc_pt_Deuteron);
  DeuteronList->Add(tof_pt_Deuteron);
  DeuteronList->Add(tpcsignal_pt_Deuteron);
  DeuteronList->Add(tpc_p_Deuteron);
  DeuteronList->Add(tof_p_Deuteron);
  DeuteronList->Add(tpcsignal_p_Deuteron);
  DeuteronList->Write("DeuteronList", TObject::kSingleKey);

  TList *ProtonList = new TList();

  ProtonList->Add(ptProton);
  ProtonList->Add(phiProton);
  ProtonList->Add(etaProton);
  ProtonList->Add(dcazProton);
  ProtonList->Add(dcaxyProton);
  ProtonList->Add(dcaz_pt_Proton);
  ProtonList->Add(dcaxy_pt_Proton);
  ProtonList->Add(tpcProton);
  ProtonList->Add(tofProton);
  ProtonList->Add(tpcsignalProton);
  ProtonList->Add(tpc_pt_Proton);
  ProtonList->Add(tof_pt_Proton);
  ProtonList->Add(tpcsignal_pt_Proton);
  ProtonList->Add(tpc_p_Proton);
  ProtonList->Add(tof_p_Proton);
  ProtonList->Add(tpcsignal_p_Proton);
  ProtonList->Write("ProtonList", TObject::kSingleKey);

  TList *LambdaList = new TList();

  LambdaList->Add(HistInvMassLambda);
  LambdaList->Add(ptLambda);
  LambdaList->Add(etaLambda);
  LambdaList->Add(phiLambda);
  LambdaList->Add(daughDCALambda);
  LambdaList->Add(transradiusLambda);
  // LambdaList->Add(HistInvMassAntiLambda);
  // LambdaList->Add(HistInvMassKaon);
  // LambdaList->Add(HistDaughterPt);
  LambdaList->Write("LambdaList", TObject::kSingleKey);

  TList *EventList = new TList();
  EventList->Add(histPosz);
  EventList->Add(histMult);
  EventList->Write("EventList", TObject::kSingleKey);

  Output->Close();
  file->Close();

  return 0;
}

Float_t ConvertBin(Char_t Input) {

  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 6.35;
  static constexpr float binned_min = -6.35;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;

  Int_t ConvInput = static_cast<Int_t>(Input);
  Float_t Output = 0.;

  if (ConvInput < underflowBin) {
    Output = binned_min;
  } else if (ConvInput > overflowBin) {
    Output = binned_max;
  } else if (ConvInput > 0) {
    Output = (ConvInput - 0.5f) * bin_width;
  } else {
    Output = (ConvInput + 0.5f) * bin_width;
  }

  return Output;
}
