// For additional jet veto
// Reference: PRL 105, 092002 (2010)
// Implementing equation number 9 from the above paper.

#ifndef NJETTINESS_H
#define NJETTINESS_H
class NJettiness {
public:

	NJettiness() {};
	~NJettiness() {};
	
    double GeneralizedTaunN(edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<edm::View<pat::Jet> > jets, double Q, double Y);
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
#endif