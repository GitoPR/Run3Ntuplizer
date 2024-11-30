//small change 
{
TCanvas *c1 =new TCanvas("c1", " ", 0, 0,700,1000);
c1->Range(0,0,1,1);
c1->SetFillColor(0);
c1->SetBorderMode(0);
c1->SetBorderSize(2);
c1->SetFrameBorderMode(0);
c1->Draw();
c1->SetLogy();
c1->SetGrid();
gStyle->SetOptStat(0);

TFile *g0 =TFile::Open("zerobias.root");
TH1F *h0 = (TH1F*)g0->Get("l1NtupleProducer/mjetcluster_pt");
TH1F *h1 = (TH1F*)g0->Get("l1NtupleProducer/nEvents");
 
int nBins = h0->GetNbinsX()+2;
 
TH1F* ratesHist0 = new TH1F("", "", nBins, h0->GetMinimum(), h0->GetMaximum());
ratesHist0->Sumw2();
int Sum0=0;

for (int i = nBins; i > 0; i--){
        Sum0 += h0->GetBinContent(i);
        ratesHist0->SetBinContent(i, Sum0);
}

float firstBin0 = ratesHist0->GetBinContent(1);
ratesHist0->Scale((double) 1.00 / firstBin0);
ratesHist0->Scale((h0->GetEntries()/h1->GetEntries()) * 40.0 * 1000000.0 / 1000.0);
ratesHist0->GetXaxis()->SetTitle("L1 jet p_{T} [GeV]");
ratesHist0->GetXaxis()->SetTitleSize(0.045);
ratesHist0->GetYaxis()->SetTitle("Rate [kHz]");
ratesHist0->GetYaxis()->SetTitleSize(0.045);
ratesHist0->SetLineWidth(2.);
ratesHist0->SetLineColor(kBlack);
ratesHist0->Draw();

TLegend *legend1 = new TLegend(0.45, 0.4, 0.7, 0.55);
legend1->SetTextFont(42);
legend1->SetLineColor(0);
legend1->SetTextSize(0.04);
legend1->SetFillColor(0);
legend1->AddEntry(ratesHist0, "cluster jet", "l");
legend1->Draw();

TLatex *t2a = new TLatex(0.5,0.9," #bf{CMS} #it{Preliminary}                     (13 TeV)");
t2a->SetNDC();
t2a->SetTextFont(42);
t2a->SetTextSize(0.045);
t2a->SetTextAlign(20);
t2a->Draw("same");
}
