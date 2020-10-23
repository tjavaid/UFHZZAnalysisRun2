git cms-init
git clone -b 10_6_12 https://ferrico@github.com/ferrico/UFHZZAnalysisRun2.git
git clone -b RunII_v2 https://github.com/VBF-HZZ/UFHZZAnalysisRun2-Accessary.git
mv UFHZZAnalysisRun2-Accessary/* ./
rm -rf UFHZZAnalysisRun2-Accessary
git clone -b tmp_Ferrico https://github.com/ferrico/KinZfitter.git
scram b -j 8
