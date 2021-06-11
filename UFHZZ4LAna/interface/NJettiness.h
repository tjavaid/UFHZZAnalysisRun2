// For additional jet veto
// Reference: PRL 105, 092002 (2010)
// Implementing equation number 9 from the above paper.

#ifndef NJETTINESS_H
#define NJETTINESS_H
class NJettiness {
public:

	NJettiness() {};
	~NJettiness() {};
	
    double GeneralizedTaunN(
        edm::Handle<pat::PackedCandidateCollection> pfcands,
        edm::Handle<edm::View<pat::Jet> > jets,
        double Q,
        double Y
        );

    void GetRapidityWeightedValues(
        unsigned int NJettiness,
        std::vector<pat::Jet> goodJets,
        float HiggsRapidity,
        float &TauB_Inc_0j,
        float &TauB_JetConstituents_0j,
        float &TauB_Inc_0j_CorrRapidity,
        float &TauB_JetConstituents_0j_CorrRapidity,
        float &TauC_Inc_0j,
        float &TauC_JetConstituents_0j,
        float &TauCnoHRapidity_Inc_0j,
        float &TauCnoHRapidity_JetConstituents_0j,
        float &TauC_Inc_0j_CorrRapidity,
        float &TauC_JetConstituents_0j_CorrRapidity,
        float &TauCnoHRapidity_Inc_0j_CorrRapidity,
        float &TauCnoHRapidity_JetConstituents_0j_CorrRapidity,
        unsigned int nJettinessSize_temp
        );
};

#endif

#ifndef NJETTINESS_CC
#define NJETTINESS_CC

/**
 * @brief      This function calculates the NJettiness
 *
 * @param[in]  pfcands  Packed candidate collection
 * @param[in]  jets     Pat jet collection
 * @param[in]  Q        square root of total invariant mass of the selected particles. Here it should be the invariant mass of 4 selected leptons.
 * @param[in]  Y        Rapidity of the vector sum of selected particles (4 selected leptons).
 */
double NJettiness::GeneralizedTaunN(edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<edm::View<pat::Jet> > jets, double Q, double Y)
{
    double tauN = -999.0;
    // FIXME: Remember to remove the selected particles from the collection of pfCandidates.
    for (const pat::PackedCandidate &pfc : *pfcands)
    {
        double dA = (pfc.pt()*TMath::Exp(Y-pfc.eta()))/Q;
        double dB = (pfc.pt()*TMath::Exp(-Y+pfc.eta()))/Q;

        tauN = TMath::Min(dA,dB);

        for (const pat::Jet &jet : *jets)
        {
            double dEta = pfc.eta() - jet.eta();
            double dPhi = pfc.phi() - jet.phi();
            double dJ = ( pfc.pt() * ( 2*TMath::CosH(dEta) - 2*TMath::Cos(dPhi) ) )/Q;
            tauN = TMath::Min(tauN, dJ);
        }   // END: for (const pat::Jet &jet : *jets)
    }   // END: for (const pat::PackedCandidate &pfc : *pfcands)

    return tauN;
}

/**
 * @brief      Gets the rapidity weighted values.
 *
 * @param[in]  NJettiness                                       N-Jettiness that we want. 0-jettiness, 1 jettiness, etc.
 *
 * @param[in]  goodJets                                         GoodJets collection
 * @param[in]  HiggsRapidity                                    The higgs rapidity
 * @param      TauB_Inc_0j                                      TauB value using inclusive jet information.
 * @param      TauB_JetConstituents_0j                          TauB value using jet constituents information.
 * @param      TauB_Inc_0j_CorrRapidity                         Rapidity corresponds to the selected jet for TauB inclusive case
 * @param      TauB_JetConstituents_0j_CorrRapidity             Rapidity corresponds to the selected jet for TauB using jet constituents
 * @param      TauC_Inc_0j                                      TauC using inclusive jet information
 * @param      TauC_JetConstituents_0j                          TauC using Jet constituents
 * @param      TauCnoHRapidity_Inc_0j                           TauC using inclusive jet, without Higgs rapidity
 * @param      TauCnoHRapidity_JetConstituents_0j               TauC using jet constituents, without Higgs rapidity
 * @param      TauC_Inc_0j_CorrRapidity                         Rapidity corresponds to selected TauC using inclusive jet information
 * @param      TauC_JetConstituents_0j_CorrRapidity             Rapidity corresponds to selected TauC using Jet constituents
 * @param      TauCnoHRapidity_Inc_0j_CorrRapidity              Rapidity corresponds to selected TauC using inclusive jet, without Higgs rapidity
 * @param      TauCnoHRapidity_JetConstituents_0j_CorrRapidity  Rapidity corresponds to selected TauC using jet constituents, without Higgs rapidity
 * @param[in]  nJettinessSize_temp                              Maximum number of jets to be considered for jettiness calculation.
 */
void NJettiness::GetRapidityWeightedValues(
        unsigned int NJettiness,
        std::vector<pat::Jet> goodJets,
        float HiggsRapidity,
        float &TauB_Inc_0j,
        float &TauB_JetConstituents_0j,
        float &TauB_Inc_0j_CorrRapidity,
        float &TauB_JetConstituents_0j_CorrRapidity,
        float &TauC_Inc_0j,
        float &TauC_JetConstituents_0j,
        float &TauCnoHRapidity_Inc_0j,
        float &TauCnoHRapidity_JetConstituents_0j,
        float &TauC_Inc_0j_CorrRapidity,
        float &TauC_JetConstituents_0j_CorrRapidity,
        float &TauCnoHRapidity_Inc_0j_CorrRapidity,
        float &TauCnoHRapidity_JetConstituents_0j_CorrRapidity,
        unsigned int nJettinessSize_temp = 0
        )
{
    unsigned int nJettinessSize;
    if (nJettinessSize_temp != 0)
    {
        if (goodJets.size() >= nJettinessSize_temp)
        {
            nJettinessSize = nJettinessSize_temp;
        }
        else
        {
            nJettinessSize = goodJets.size();
        }
    }
    else
    {
            nJettinessSize = goodJets.size();
    }

    float TauB_Inc_0j_temp = -999.0;
    float TauC_Inc_0j_temp = -999.0;
    float TauB_JetConstituents_0j_temp = -999.0;
    float TauC_JetConstituents_0j_temp = -999.0;
    float TauCnoHRapidity_JetConstituents_0j_temp = -999.0;

    float TauCnoHRapidity_Inc_0j_temp = -999.0;

    float TauC_Inc_0j_CorrRapidity_temp = -999.0;
    float TauB_Inc_0j_CorrRapidity_temp = -999.0;
    float TauB_JetConstituents_0j_CorrRapidity_temp = -999.0;
    float TauC_JetConstituents_0j_CorrRapidity_temp = -999.0;
    float TauCnoHRapidity_JetConstituents_0j_CorrRapidity_temp = -999.0;
    float TauCnoHRapidity_Inc_0j_CorrRapidity_temp = -999.0;

    for (unsigned int JetCounter = NJettiness; JetCounter < nJettinessSize; ++JetCounter)
    {
        // Inclusive tauC
        float TauC_Inc_0j_num = sqrt(goodJets[JetCounter].energy()*goodJets[JetCounter].energy() - goodJets[JetCounter].pz()*goodJets[JetCounter].pz());
        float TauC_Inc_0j_den = 2*cosh(goodJets[JetCounter].rapidity() - HiggsRapidity);
        if (TauC_Inc_0j_num/TauC_Inc_0j_den > TauC_Inc_0j_temp)
        {
            TauC_Inc_0j_temp = TauC_Inc_0j_num/TauC_Inc_0j_den;
            TauC_Inc_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }

        // Inclusive tauC without Higgs rapidity
        if (TauC_Inc_0j_num/(2*cosh(goodJets[JetCounter].rapidity())) > TauCnoHRapidity_Inc_0j_temp)
        {
            TauCnoHRapidity_Inc_0j_temp = TauC_Inc_0j_num/(2*cosh(goodJets[JetCounter].rapidity()));
            TauCnoHRapidity_Inc_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }

        // Inclusive tauB
        if (abs(goodJets[JetCounter].energy() - goodJets[JetCounter].pz()) > TauB_Inc_0j_temp)
        {
            TauB_Inc_0j_temp = abs(goodJets[JetCounter].energy() - goodJets[JetCounter].pz());
            TauB_Inc_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }

        // TauB & TauC using Jet Constituents
        float TauB_JetConstituents_0j_local = 0.0;
        float TauC_JetConstituents_0j_local = 0.0;
        float TauC_JetConstituents_0j_local2 = 0.0;
        for ( auto const & constituent : goodJets[JetCounter].daughterPtrVector())
        {
            // tauB
            TauB_JetConstituents_0j_local += abs(constituent->energy() - constituent->pz());

            // tauC
            double TauC2_numerator = sqrt(constituent->energy()*constituent->energy()-constituent->pz()*constituent->pz());
            // tauC: version1: with higgs rapidity
            double TauC2_denominator = 2*cosh(constituent->rapidity() - HiggsRapidity);
            TauC_JetConstituents_0j_local += TauC2_numerator/TauC2_denominator;

            // tauC: version1: without higgs rapidity
            TauC_JetConstituents_0j_local2 += TauC2_numerator/(2*cosh(constituent->rapidity()));
        }
        if (TauB_JetConstituents_0j_local > TauB_JetConstituents_0j_temp)
        {
            TauB_JetConstituents_0j_temp = TauB_JetConstituents_0j_local;
            TauB_JetConstituents_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }
        if (TauC_JetConstituents_0j_local > TauC_JetConstituents_0j_temp)
        {
            TauC_JetConstituents_0j_temp = TauC_JetConstituents_0j_local;
            TauC_JetConstituents_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }
        if (TauC_JetConstituents_0j_local2 > TauCnoHRapidity_JetConstituents_0j_temp)
        {
            TauCnoHRapidity_JetConstituents_0j_temp = TauC_JetConstituents_0j_local2;
            TauCnoHRapidity_JetConstituents_0j_CorrRapidity_temp = goodJets[JetCounter].rapidity();
        }
    }
    TauB_Inc_0j = TauB_Inc_0j_temp;
    TauB_JetConstituents_0j = TauB_JetConstituents_0j_temp;

    TauB_Inc_0j_CorrRapidity = TauB_Inc_0j_CorrRapidity_temp;
    TauB_JetConstituents_0j_CorrRapidity = TauB_JetConstituents_0j_CorrRapidity_temp;

    TauC_Inc_0j = TauC_Inc_0j_temp;
    TauC_JetConstituents_0j = TauC_JetConstituents_0j_temp;
    TauCnoHRapidity_Inc_0j = TauCnoHRapidity_Inc_0j_temp;
    TauCnoHRapidity_JetConstituents_0j = TauCnoHRapidity_JetConstituents_0j_temp;

    TauC_Inc_0j_CorrRapidity = TauC_Inc_0j_CorrRapidity_temp;
    TauC_JetConstituents_0j_CorrRapidity = TauC_JetConstituents_0j_CorrRapidity_temp;
    TauCnoHRapidity_Inc_0j_CorrRapidity = TauCnoHRapidity_Inc_0j_CorrRapidity_temp;
    TauCnoHRapidity_JetConstituents_0j_CorrRapidity = TauCnoHRapidity_JetConstituents_0j_CorrRapidity_temp;
}
#endif