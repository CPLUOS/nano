cd $CMSSW_BASE/src/nano/
wget http://www.hep.phy.cam.ac.uk/~lester/dtm662/mt2/Releases/oxbridgekinetics.tar.gz
tar -xzvf oxbridgekinetics.tar.gz
rm oxbridgekinetics.tar.gz
mv oxbridgekinetics-1.2/ oxbridgekinetics
cd oxbridgekinetics
touch BuildFile.xml
#####################################################################
#make a buildfile for integrating oxbridgekinetics library into nano
#specify flags in BuildFile to ignore unnecessary warnings and errors
cat > BuildFile.xml << EOL
<use name="root" />
<flags LDFLAGS="-lMinuit2" />
<flags CXXFLAGS="-Wno-error=unused-but-set-variable" />
<flags CXXFLAGS="-Wno-error=unused-variable" />
<flags CXXFLAGS="-Wno-error=sign-compare" />
<flags CXXFLAGS="-Wno-error=maybe-uninitialized" />
<export>
  <lib name="1"/>
</export>
EOL
#####################################################################
mv Mt2/ src
cd src/Mt2
sed -i s:"Mt2/:": * # modify the path of Mt2 header files.
cd ../../
scram b
#set right path in the nano/analysis/bin/BuildFile.xml
sed -i '/erase/d' $CMSSW_BASE/src/nano/analysis/bin/BuildFile.xml
