#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# File              : postProcessing.py
# Author            : Anton Riedel <anton.riedel@tum.de>
# Date              : 10.11.2022
# Last Modified Date: 02.12.2022
# Last Modified By  : Anton Riedel <anton.riedel@tum.de>

import ROOT
import json
import numpy as np
import uproot
import sys

Particles = [
    "Proton",
    "Deuteron",
    "Lambda",
    "PosDaughter",
    "NegDaughter",
]

RawParticles = ["RawTrack", "RawLambda", "RawPosDaughter", "RawNegDaughter"]


def SetupHistograms(HistConfigFileName):
    # load config for all histograms
    Hists = {}
    TempHist = {}
    with open(HistConfigFileName, "r") as HistConfigFile:
        HistConfig = json.load(HistConfigFile)
        for key in HistConfig["Event"]:
            TempHist[key] = ROOT.TH1F(
                "Event_" + key,
                "Event_" + key,
                HistConfig["Event"][key]["Bins"],
                HistConfig["Event"][key]["RangeMin"],
                HistConfig["Event"][key]["RangeMax"],
            )
        Hists["Event"] = TempHist
        for P in Particles + RawParticles:
            TempHist = {}
            for key in HistConfig["Particle_1D"]:
                TempHist[key] = ROOT.TH1F(
                    P + "_" + key,
                    P + "_" + key,
                    HistConfig["Particle_1D"][key]["Bins"],
                    HistConfig["Particle_1D"][key]["RangeMin"],
                    HistConfig["Particle_1D"][key]["RangeMax"],
                )
            for key in HistConfig["Particle_2D"]:
                TempHist[key] = ROOT.TH2F(
                    P + "_" + key,
                    P + "_" + key,
                    HistConfig["Particle_2D"][key]["XBins"],
                    HistConfig["Particle_2D"][key]["XRangeMin"],
                    HistConfig["Particle_2D"][key]["XRangeMax"],
                    HistConfig["Particle_2D"][key]["YBins"],
                    HistConfig["Particle_2D"][key]["YRangeMin"],
                    HistConfig["Particle_2D"][key]["YRangeMax"],
                )
            Hists[P] = TempHist
    return Hists


def SetupCuts(ConfigFileName):
    # load config for all histograms
    Cuts = {}
    TempCuts = {}
    with open(ConfigFileName, "r") as ConfigFile:
        Config = json.load(ConfigFile)
        for P in Particles + ["Event"]:
            TempCuts = {}
            for key in Config[P]:
                TempCuts[key] = (
                    Config[P][key]["Min"],
                    Config[P][key]["Max"],
                )
                if P == "Proton":
                    TempCuts["Proton_PTPC"] = Config["Proton_PTPC"]
            Cuts[P] = TempCuts
    return Cuts


def CheckCut(Cut, CutKey, Tree, TreeKey, TreeIndex):
    if Cut[CutKey][0] <= Tree[TreeKey][TreeIndex] <= Cut[CutKey][1]:
        return True
    else:
        return False


def CheckEvent(CollisionCuts, CollisionTree, Index):
    if CheckCut(CollisionCuts, "VertexZ", CollisionTree, "fPosZ", Index):
        return True
    else:
        return False


def CheckProton(ProtonCuts, ParticleTree, ParticleDebugTree, Index):
    TPCCrossedRowsOverFindable = float(
        ParticleDebugTree["fTPCNClsCrossedRows"][Index]
    ) / float(ParticleDebugTree["fTPCNClsFindable"][Index])

    if (
        CheckCut(ProtonCuts, "Charge", ParticleDebugTree, "fSign", Index)
        and CheckCut(ProtonCuts, "Pt", ParticleTree, "fPt", Index)
        and CheckCut(ProtonCuts, "Eta", ParticleTree, "fEta", Index)
        and CheckCut(ProtonCuts, "DCAz", ParticleDebugTree, "fDcaZ", Index)
        and CheckCut(ProtonCuts, "DCAxy", ParticleDebugTree, "fDcaXY", Index)
        and CheckCut(
            ProtonCuts, "TPCClustersFound", ParticleDebugTree, "fTPCNClsFound", Index
        )
        and CheckCut(
            ProtonCuts,
            "TPCCrossedRows",
            ParticleDebugTree,
            "fTPCNClsCrossedRows",
            Index,
        )
        and ProtonCuts["TPCCrossedRowsOverFindable"][0]
        <= TPCCrossedRowsOverFindable
        <= ProtonCuts["TPCCrossedRowsOverFindable"][1]
        and CheckCut(
            ProtonCuts, "TPCClustersShared", ParticleDebugTree, "fTPCNClsShared", Index
        )
    ):
        P = ParticleTree["fPt"][Index] * np.cosh(ParticleTree["fEta"][Index])
        TPC = ConvertBin(ParticleDebugTree["fTPCNSigmaStorePr"][Index])
        TOF = ConvertBin(ParticleDebugTree["fTOFNSigmaStorePr"][Index])
        if P <= ProtonCuts["Proton_PTPC"]:
            if ProtonCuts["NSigmaTPC"][0] <= TPC <= ProtonCuts["NSigmaTPC"][1]:
                return True
            else:
                return False
        else:
            if (
                ProtonCuts["NSigmaTPCTOF"][0]
                <= np.sqrt(TPC**2 + TOF**2)
                <= ProtonCuts["NSigmaTPCTOF"][1]
            ):
                return True
            else:
                return False
    else:
        return False


def CheckDeuteron(DeuteronCuts, ParticleTree, ParticleDebugTree, Index):
    TPCCrossedRowsOverFindable = float(
        ParticleDebugTree["fTPCNClsCrossedRows"][Index]
    ) / float(ParticleDebugTree["fTPCNClsFindable"][Index])
    TPC = ConvertBin(ParticleDebugTree["fTPCNSigmaStoreDe"][Index])

    if (
        CheckCut(DeuteronCuts, "Charge", ParticleDebugTree, "fSign", Index)
        and CheckCut(DeuteronCuts, "Pt", ParticleTree, "fPt", Index)
        and CheckCut(DeuteronCuts, "Eta", ParticleTree, "fEta", Index)
        and CheckCut(DeuteronCuts, "DCAz", ParticleDebugTree, "fDcaZ", Index)
        and CheckCut(DeuteronCuts, "DCAxy", ParticleDebugTree, "fDcaXY", Index)
        and CheckCut(
            DeuteronCuts, "TPCClustersFound", ParticleDebugTree, "fTPCNClsFound", Index
        )
        and CheckCut(
            DeuteronCuts,
            "TPCCrossedRows",
            ParticleDebugTree,
            "fTPCNClsCrossedRows",
            Index,
        )
        and DeuteronCuts["TPCCrossedRowsOverFindable"][0]
        <= TPCCrossedRowsOverFindable
        <= DeuteronCuts["TPCCrossedRowsOverFindable"][1]
        and CheckCut(
            DeuteronCuts,
            "TPCClustersShared",
            ParticleDebugTree,
            "fTPCNClsShared",
            Index,
        )
        and CheckCut(DeuteronCuts, "ITSClusters", ParticleDebugTree, "fITSNCls", Index)
        and CheckCut(
            DeuteronCuts,
            "ITSClustersIB",
            ParticleDebugTree,
            "fITSNClsInnerBarrel",
            Index,
        )
        and DeuteronCuts["NSigmaTPC"][0] <= TPC <= DeuteronCuts["NSigmaTPC"][1]
    ):
        TPC_Proton = ConvertBin(ParticleDebugTree["fTPCNSigmaStorePr"][Index])
        TPC_Pion = ConvertBin(ParticleDebugTree["fTPCNSigmaStorePi"][Index])
        TPC_Electron = ConvertBin(ParticleDebugTree["fTPCNSigmaStoreEl"][Index])
        if (
            DeuteronCuts["TPCRejection"][0]
            <= TPC_Proton
            <= DeuteronCuts["TPCRejection"][1]
            or DeuteronCuts["TPCRejection"][0]
            <= TPC_Pion
            <= DeuteronCuts["TPCRejection"][1]
            or DeuteronCuts["TPCRejection"][0]
            <= TPC_Electron
            <= DeuteronCuts["TPCRejection"][1]
        ):
            return False
        else:
            return True
    else:
        return False


def CheckLambda(LambdaCuts, ParticleTree, ParticleDebugTree, Index, PrimaryVertex):
    DecayVertexDist = np.sqrt(
        (ParticleDebugTree["fDecayVtxX"][Index] - PrimaryVertex[0]) ** 2
        + (ParticleDebugTree["fDecayVtxY"][Index] - PrimaryVertex[1]) ** 2
        + (ParticleDebugTree["fDecayVtxZ"][Index] - PrimaryVertex[2]) ** 2
    )
    if (
        CheckCut(LambdaCuts, "Pt", ParticleTree, "fPt", Index)
        and CheckCut(LambdaCuts, "Eta", ParticleTree, "fEta", Index)
        and CheckCut(LambdaCuts, "CosPA", ParticleTree, "fTempFitVar", Index)
        and CheckCut(
            LambdaCuts, "TransRadius", ParticleDebugTree, "fTransRadius", Index
        )
        and LambdaCuts["DecayVertexDist"][0]
        <= DecayVertexDist
        <= LambdaCuts["DecayVertexDist"][1]
        and CheckCut(LambdaCuts, "DaughterDCA", ParticleDebugTree, "fDaughDCA", Index)
        and CheckCut(LambdaCuts, "LambdaInvMass", ParticleTree, "fMLambda", Index)
        and not CheckCut(LambdaCuts, "K0InvMass", ParticleDebugTree, "fMKaon", Index)
    ):
        return True
    else:
        return False


def CheckDaugher(DaughterCuts, ParticleTree, ParticleDebugTree, Index, TPCNSigmaBranch):
    DCAPrimaryVertex = np.sqrt(
        ParticleDebugTree["fDcaXY"][Index] ** 2 + ParticleDebugTree["fDcaZ"][Index] ** 2
    )
    TPC = ConvertBin(ParticleDebugTree[TPCNSigmaBranch][Index])
    if (
        CheckCut(DaughterCuts, "Charge", ParticleDebugTree, "fSign", Index)
        and CheckCut(DaughterCuts, "Eta", ParticleTree, "fEta", Index)
        and CheckCut(
            DaughterCuts, "TPCClustersFound", ParticleDebugTree, "fTPCNClsFound", Index
        )
        and DaughterCuts["NSigmaTPC"][0] <= TPC <= DaughterCuts["NSigmaTPC"][1]
        and not (
            DaughterCuts["DCAPrimaryVertex"][0]
            <= DCAPrimaryVertex
            <= DaughterCuts["DCAPrimaryVertex"][1]
        )
    ):
        return True
    else:
        return False


def ConvertBin(Input):
    Size = 1  # sizeof(int8_t)
    nbins = (1 << 8 * Size) - 2
    overflowBin = nbins >> 1
    underflowBin = -(nbins >> 1)
    binned_max = 6.35
    binned_min = -6.35
    bin_width = (binned_max - binned_min) / nbins
    Output = 0.0
    if Input < underflowBin:
        Output = binned_min
    elif Input > overflowBin:
        Output = binned_max
    elif Input > 0:
        Output = (Input - 0.5) * bin_width
    else:
        Output = (Input + 0.5) * bin_width
    return Output


def ProcessTrack(
    TreeIndex,
    ParticleTree,
    ParticleDebugTree,
    Histograms,
    TPCNSigmaBranch,
    TOFNSigmaBranch,
    PrimaryVertex,
):

    # compute quantities that cannot be pulled from tree
    if ParticleDebugTree["fTPCNClsFindable"][TreeIndex] != 0:
        TPCCrossedRowsOverFindable = float(
            ParticleDebugTree["fTPCNClsCrossedRows"][TreeIndex]
        ) / float(ParticleDebugTree["fTPCNClsFindable"][TreeIndex])
    else:
        TPCCrossedRowsOverFindable = 3
    P = ParticleTree["fPt"][TreeIndex] * np.cosh(ParticleTree["fEta"][TreeIndex])
    if TPCNSigmaBranch != "":
        TPC = ConvertBin(ParticleDebugTree[TPCNSigmaBranch][TreeIndex])
        TOF = ConvertBin(ParticleDebugTree[TOFNSigmaBranch][TreeIndex])
    else:
        TPC = 0
        TOF = 0
    DCAPrimaryVertex = np.sqrt(
        ParticleDebugTree["fDcaXY"][TreeIndex] ** 2
        + ParticleDebugTree["fDcaZ"][TreeIndex] ** 2
    )
    DecayVertexDist = np.sqrt(
        (ParticleDebugTree["fDecayVtxX"][TreeIndex] - PrimaryVertex[0]) ** 2
        + (ParticleDebugTree["fDecayVtxY"][TreeIndex] - PrimaryVertex[1]) ** 2
        + (ParticleDebugTree["fDecayVtxZ"][TreeIndex] - PrimaryVertex[2]) ** 2
    )

    # fill 1D histograms
    Histograms["Charge"].Fill(ParticleDebugTree["fSign"][TreeIndex])
    Histograms["Pt"].Fill(ParticleTree["fPt"][TreeIndex])
    Histograms["Eta"].Fill(ParticleTree["fEta"][TreeIndex])
    Histograms["Phi"].Fill(ParticleTree["fPhi"][TreeIndex])
    Histograms["DCAxy"].Fill(ParticleDebugTree["fDcaXY"][TreeIndex])
    Histograms["DCAz"].Fill(ParticleDebugTree["fDcaZ"][TreeIndex])
    Histograms["DCAPrimaryVertex"].Fill(DCAPrimaryVertex)
    Histograms["DaughterDCA"].Fill(ParticleDebugTree["fDaughDCA"][TreeIndex])
    Histograms["TPCClustersFound"].Fill(ParticleDebugTree["fTPCNClsFound"][TreeIndex])
    Histograms["TPCClustersFindable"].Fill(
        ParticleDebugTree["fTPCNClsFindable"][TreeIndex]
    )
    Histograms["TPCCrossedRows"].Fill(
        ParticleDebugTree["fTPCNClsCrossedRows"][TreeIndex]
    )
    Histograms["TPCCrossedRowsOverFindable"].Fill(TPCCrossedRowsOverFindable)
    Histograms["TPCClustersShared"].Fill(ParticleDebugTree["fTPCNClsShared"][TreeIndex])
    Histograms["ITSClusters"].Fill(ParticleDebugTree["fITSNCls"][TreeIndex])
    Histograms["ITSClustersIB"].Fill(
        ParticleDebugTree["fITSNClsInnerBarrel"][TreeIndex]
    )
    Histograms["NSigmaTPC"].Fill(TPC)
    Histograms["NSigmaTOF"].Fill(TOF)
    Histograms["CosPA"].Fill(ParticleTree["fTempFitVar"][TreeIndex])
    Histograms["TransRadius"].Fill(ParticleDebugTree["fTransRadius"][TreeIndex])
    Histograms["DecayVertexDist"].Fill(DecayVertexDist)
    Histograms["K0InvMass"].Fill(ParticleDebugTree["fMKaon"][TreeIndex])
    Histograms["LambdaInvMass"].Fill(ParticleTree["fMLambda"][TreeIndex])

    # fill 2D Histograms
    Histograms["DCAzVsPt"].Fill(
        ParticleTree["fPt"][TreeIndex], ParticleDebugTree["fDcaZ"][TreeIndex]
    )
    Histograms["DCAxyVsPt"].Fill(
        ParticleTree["fPt"][TreeIndex], ParticleDebugTree["fDcaXY"][TreeIndex]
    )
    Histograms["NSigmaTPCvsP"].Fill(P, TPC)
    Histograms["NSigmaTOFvsP"].Fill(P, TOF)


def main(InputFileName, OutputFile, HistConfigFileName, ConfigFileName):

    # setup histograms
    Hists = SetupHistograms(HistConfigFileName)
    Cuts = SetupCuts(ConfigFileName)

    ProcessedCollisions = set()
    CollisionIndex = -1

    # open root file with uproot
    with uproot.open(InputFileName) as file:

        # get all top level TDirectoryFile
        Directories = [
            Dir for Dir in file.classnames() if not "/" in Dir and Dir.startswith("DF_")
        ]

        # now loop through the top level TDirectoryFiles
        for Dir in Directories:

            # and get the trees
            TreeParticle = file[Dir + "/O2femtodreamparts;1"].arrays(library="np")
            TreeParticleDebug = file[Dir + "/O2femtodebugparts;1"].arrays(library="np")
            TreeEvents = file[Dir + "/O2femtodreamcols;1"].arrays(library="np")

            Entries = np.shape(TreeParticle["fPartType"])[0]

            # loop through the trees
            for TreeIndex in range(Entries):

                CollisionIndex = TreeParticle["fIndexFemtoDreamCollisions"][TreeIndex]
                PrimaryVertex = [0, 0, TreeEvents["fPosZ"][CollisionIndex]]

                # cut event
                if not CheckEvent(Cuts["Event"], TreeEvents, CollisionIndex):
                    continue

                if CollisionIndex not in ProcessedCollisions:
                    Hists["Event"]["VertexZ"].Fill(TreeEvents["fPosZ"][CollisionIndex])
                    Hists["Event"]["Multiplicity"].Fill(
                        TreeEvents["fMultV0M"][CollisionIndex]
                    )
                    ProcessedCollisions.add(CollisionIndex)

                # if the particle is a track, fill proton and deuteron histsograms
                if TreeParticle["fPartType"][TreeIndex] == 0:
                    ProcessTrack(
                        TreeIndex,
                        TreeParticle,
                        TreeParticleDebug,
                        Hists["RawTrack"],
                        "",
                        "",
                        PrimaryVertex,
                    )
                    if CheckProton(
                        Cuts["Proton"], TreeParticle, TreeParticleDebug, TreeIndex
                    ):
                        ProcessTrack(
                            TreeIndex,
                            TreeParticle,
                            TreeParticleDebug,
                            Hists["Proton"],
                            "fTPCNSigmaStorePr",
                            "fTOFNSigmaStorePr",
                            PrimaryVertex,
                        )
                    if CheckDeuteron(
                        Cuts["Deuteron"], TreeParticle, TreeParticleDebug, TreeIndex
                    ):
                        ProcessTrack(
                            TreeIndex,
                            TreeParticle,
                            TreeParticleDebug,
                            Hists["Deuteron"],
                            "fTPCNSigmaStoreDe",
                            "fTOFNSigmaStoreDe",
                            PrimaryVertex,
                        )
                # if the particle is a v0, fill lambda histograms
                elif TreeParticle["fPartType"][TreeIndex] == 1:
                    # protect against running out of bounds
                    if TreeIndex + 2 >= Entries:
                        break

                    ProcessTrack(
                        TreeIndex,
                        TreeParticle,
                        TreeParticleDebug,
                        Hists["RawLambda"],
                        "",
                        "",
                        PrimaryVertex,
                    )
                    ProcessTrack(
                        TreeIndex + 1,
                        TreeParticle,
                        TreeParticleDebug,
                        Hists["RawPosDaughter"],
                        "fTPCNSigmaStorePr",
                        "fTOFNSigmaStorePr",
                        PrimaryVertex,
                    )
                    ProcessTrack(
                        TreeIndex + 2,
                        TreeParticle,
                        TreeParticleDebug,
                        Hists["RawNegDaughter"],
                        "fTPCNSigmaStorePi",
                        "fTOFNSigmaStorePi",
                        PrimaryVertex,
                    )
                    if (
                        CheckLambda(
                            Cuts["Lambda"],
                            TreeParticle,
                            TreeParticleDebug,
                            TreeIndex,
                            PrimaryVertex,
                        )
                        and CheckDaugher(
                            Cuts["PosDaughter"],
                            TreeParticle,
                            TreeParticleDebug,
                            TreeIndex + 1,
                            "fTPCNSigmaStorePr",
                        )
                        and CheckDaugher(
                            Cuts["NegDaughter"],
                            TreeParticle,
                            TreeParticleDebug,
                            TreeIndex + 2,
                            "fTPCNSigmaStorePi",
                        )
                    ):
                        ProcessTrack(
                            TreeIndex,
                            TreeParticle,
                            TreeParticleDebug,
                            Hists["Lambda"],
                            "",
                            "",
                            PrimaryVertex,
                        )
                        ProcessTrack(
                            TreeIndex + 1,
                            TreeParticle,
                            TreeParticleDebug,
                            Hists["PosDaughter"],
                            "fTPCNSigmaStorePr",
                            "fTOFNSigmaStorePr",
                            PrimaryVertex,
                        )
                        ProcessTrack(
                            TreeIndex + 2,
                            TreeParticle,
                            TreeParticleDebug,
                            Hists["NegDaughter"],
                            "fTPCNSigmaStorePi",
                            "fTOFNSigmaStorePi",
                            PrimaryVertex,
                        )

    # save result into root file
    OutputFile = ROOT.TFile(OutputFileName, "recreate")
    for P in Particles + RawParticles + ["Event"]:
        List = ROOT.TList()
        for hist in Hists[P].values():
            List.Add(hist)
        List.Write(P, 1)  # 1 = TObject::kSingleKey
    OutputFile.Close()


if __name__ == "__main__":

    # input handling
    InputFileName = sys.argv[1]
    OutputFileName = sys.argv[2]
    HistConfigFileName = sys.argv[3]
    ConfigFileName = sys.argv[4]

    # hardcode values for testing
    # InputFileName = "../FemtoAO2D.root"
    # OutputFileName = "test.root"
    # HistConfigFileName = "./HistConfig.json"
    # ConfigFileName = "./StandardCuts.json"

    main(InputFileName, OutputFileName, HistConfigFileName, ConfigFileName)
