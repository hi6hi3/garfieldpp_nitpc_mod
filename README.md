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
$ make
$ cd ..
$ ./nitpc
```
