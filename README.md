---
# garfieldpp_nitpc
## 陰イオンガスTPCのためのGarfield++シミュレーションツール
---
## source code author 
Takuya Shimada

---
## How to
.bashrc内で環境変数を定義しておく
```
# .bashrc内
export GARFIELD_HOME=/pash/to/garfieldpp_nitpc
export HEED_DATABASE=$GARFIELD_HOME/Heed/heed++/databece
```

```
$ cd $GARFIELD_HOME
$ make
```
---
## Example
```
$ cd $GARFIELD_HOME/Example/NITPC/src
$ vim nitpc.cxx
```
nitpc.cxxを編集
```
// nitpc.cxx内
...
//gas->LoadIonMobility("/work/shimada/Garfield++/Mobility/IonMobility_SF6-_SF6.txt"); コメントアウト
gas->LoadIonMobility("/path/to/Data/IonMobility_SF6-_SF6.txt"); // Mobilityのデータの絶対パスを指定
...
//std::string data_dir = "/work/shimada/Garfield++/GEM"; コメントアウト
std::string data_dir = "/path/to/gem/"; // Elmerからのアウトプットのディレクトリにパスを通す
std::string header = data_dir + "/gemcell/mesh.header";
std::string element = data_dir + "/gemcell/mesh.elements";
std::string node = data_dir + "/gemcell/mesh.nodes";
std::string eps = data_dir + "/gemcell/dielectrics.dat";
std::string volt = data_dir + "/gemcell/gemcell.result";
...
```
コンパイルして実効すると、シミュレーションが走る
```
$ make
$ cd ..
$ ./nitpc
....
create track.png
....
$ display track.png
```
track.pngを見ると、陰イオン(青)と電子(橙)がドリフト・ガス増幅している事がわかる

時間経過を見る
![nitpc gem simulation](https://ppwww.phys.sci.kobe-u.ac.jp/~newage/Shimada/share/avalanche_negativeion.gif)

