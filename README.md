***
# CFD _paricle _method _basics
***
<br>

粒子法を用いた簡単な水柱崩壊のプログラム（非圧縮流体）

# このプログラムについて
このプログラムでは、粒子法を用いた水柱崩壊のシミュレーションを行う。以下の参考書の付録のc++で書かれたプログラムを、FORTRANで使えるようにすることを主な目的としている。また、この文章は自分の頭の中を整理するために書いていると言う節もある。

## 参考書
粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜,丸善出版

<br><br><br>
***
# 粒子法の概要
***
<br>

以下は、「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」の内容の要点を抜き出したものである。基本的には、微分方程式を時間ステップについて離散化して各時間ステップごとに粒子の位置を更新and出力していくという、通常の微分方程式の数値解法の考え方と大差はない。ただし空間gridを切らないという点で通常の格子法とは空間の離散化方法が大きく異なる。

## Navier-stokes equation
ラグランジュ形式で書いたNavier-stokes equationは以下のよう

$$
\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}
$$

ただし $\nu$ は一定とし、また非圧縮流体を考え $\rho = \rho^0 = \text{const}$ と考えることにする。以下ではこの方程式に従う流体の運動を数値的に解くことを考える。

## 粒子法という考え方
粒子法においては、流体を微小な流体粒子の集まりとして表現し、上記のナビエストークス方程式を各粒子が存在する位置において解くことで粒子の位置を更新する。すなわち上の式を各粒子 $i$ の位置で解き、

$$
\begin{aligned}
&\boldsymbol{u}^{k+1} _i = \boldsymbol{u}^k _i + \Delta t\cdot\left(-\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g} \right)^{k} _i \newline 
&\boldsymbol{r}^{k+1} _i= \boldsymbol{r}^k _i + \Delta t\cdot  \boldsymbol{u}^{k}
\end{aligned}
$$

によって粒子の時間ステップ $k+1$ における位置を更新する。この時鍵になるのが影響半径という考え方である。すなわち、注目している粒子から（恣意的に定めて）一定の距離内のみに存在する粒子について計算を行う。より具体的には、粒子から離れるほど重みが $0$ に近づくような重み関数を計算に用いることになる。

## 微分演算子の差分化

粒子法は格子法のように単純なグリッドで空間を離散化するわけではないので、単純に差分を用いて微分演算子を離散化することはできない。そこでまずは粒子法におけるナブラ演算子（一階微分）及びラプラシアン演算子（二回微分）の離散化を考える必要がある。

### ナブラ演算子（圧力項）の離散化

まずはじめに空間方向の離散化ナブラ演算子 $\nabla$ の離散化を考える。

まずは、すごいざっくりした場合で考察する。今、二次元系で粒子がブワーってあったとする。さらに、粒子間のベクトル $\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}$ と $\boldsymbol{X} _{i+2}-\boldsymbol{X} _{i}$ が正規直行系をなすとする。

```txt
        o <- i+2
        |
   i -> o----o <- i+1
```

今ある座標系 $X,Y$ に対し一般に

$$
\nabla\phi = \dfrac{\partial \phi}{\partial X}\nabla X + \dfrac{\partial \phi}{\partial Y}\nabla Y
$$

で、 $X,Y$ が正規直交なら $\nabla X=\hat{\boldsymbol{e}} _X$ などから

$$
\nabla\phi = \dfrac{\partial \phi}{\partial X}\hat{\boldsymbol{e}} _X + \dfrac{\partial \phi}{\partial Y}\hat{\boldsymbol{e}} _Y
$$

となる。よって、

$$
\nabla\phi = \dfrac{\partial \phi}{\partial X _{\text{i+１方向}}}\dfrac{\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}|} + \dfrac{\partial \phi}{\partial X _{\text{i+2方向}}}\dfrac{\boldsymbol{X} _{i+2}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{i+2}-\boldsymbol{X} _{i}|}
$$

であり、偏微分演算子を離散化すれば

$$
\nabla\phi = \dfrac{\phi(\boldsymbol{X} _{i+1})-\phi(\boldsymbol{X} _{i})}{|\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}|}\dfrac{\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}|} + \dfrac{\phi(\boldsymbol{X} _{i+2})-\phi(\boldsymbol{X} _{i})}{|\boldsymbol{X} _{i+2}-\boldsymbol{X} _{i}|}\dfrac{\boldsymbol{X} _{i+2}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{i+1}-\boldsymbol{X} _{i}|}
$$

となる。よって今、影響半径内の総粒子数を $N$ として、全ての粒子 $\boldsymbol{X} _{j}$ に対し $\boldsymbol{X} _{j}-\boldsymbol{X} _{i}$ と $\boldsymbol{X} _{j'}-\boldsymbol{X} _{i}$　を満たすような $\boldsymbol{X} _{j'}$ が存在している（すなわち上の式が成立）と仮定する。このとき、 $N$ この粒子に対し上式を足し合わせれば（どの粒子も二回現れることに注意して）

$$
N\nabla\phi = 2\displaystyle\sum _{i\neq j}\dfrac{\phi(\boldsymbol{X} _{j})-\phi(\boldsymbol{X} _{i})}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}\dfrac{\boldsymbol{X} _{j}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}
$$

となる。

さて、実際の液体中の粒子はこんなに都合が良くはない。そもそも、綺麗に正規直行関係にある粒子ペアが必ず見つかるとは限らない上に、粒子 $X _{i}$ から遠い粒子が式に食い込むほど近似精度は悪くなってしまう。そこで、「まあざっくり上のような足し合わせでええやろ、でも重み関数をかけながら足し合わせて行こう」みたいなことを考えるわけで、それにより

$$
\left(\displaystyle\sum _{i\neq j}{w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)}\right)\nabla\phi = 2\displaystyle\sum _{i\neq j}\dfrac{\phi(\boldsymbol{X} _{j})-\phi(\boldsymbol{X} _{i})}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}\dfrac{\boldsymbol{X} _{j}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

という離散化の式を考えることができる。ここで、重み関数とは距離が遠くなるほど値が小さくなり適当な影響半径 $r _e$ で０となるような関数であり、MPS法では次の重み関数を用いる。

$$
w \left( r\right) = \left\lbrace \begin{aligned}&\dfrac{r _e}{r}-1&\quad\left(r<r _e\right)\newline
&0 &\quad\left(r>r _e\right)\end{aligned} \right.
$$

また、今、上の式の $2$ は二次元を仮定したことにより出てきた式であるから、次元数 $d$ をもちいて $2\rightarrow d $を書き換えることにする。また今、重みつき和の規格化定数である

$$
n^i = \displaystyle\sum _{i\neq j}{w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)}
$$

については、初期状態で適当な液体内部粒子について計算しておいた値

$$
n _0 = \displaystyle\sum _{i\neq j}{w(|\boldsymbol{X} _{j}^{0}-\boldsymbol{X} _{i}^{0}|)}
$$

を用いて計算することにする。これを用いれば、

$$
\nabla\phi = \dfrac{d}{n^0}\displaystyle\sum _{i\neq j}\dfrac{\phi(\boldsymbol{X} _{j})-\phi(\boldsymbol{X} _{i})}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}\dfrac{\boldsymbol{X} _{j}-\boldsymbol{X} _{i}}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|}w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

として粒子法における空間一階微分の離散化が得られる。

### ラプラシアン（粘性項）の離散化
ナブラの場合と似たように考えていくことができて（めんどくさいので割愛、詳しくは「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」を参照。）、

$$
\nabla^2 \phi= \dfrac{2d}{\lambda _0 n^0}\displaystyle\sum _{i\neq j}(\phi(\boldsymbol{X} _{j})-\phi(\boldsymbol{X} _{i}))w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

のように離散化することができる。ただし、

$$
\lambda _0 = \dfrac{1}{n^0}\displaystyle\sum _{i\neq j}|\boldsymbol{X} _j -\boldsymbol{X} _i|^2 w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

である。

## 半陰解法

微分演算子の離散化がもとまったので、実際に粒子速度・粒子位置をナビエストークス方程式に従って更新していくことを考える。今改めて、ナビエストークス方程式は以下のようであった。

$$
\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho^0}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}
$$

（ただし非圧縮条件を考え $\rho = \rho^0$ としている。）

さて、右辺第二項及び第三項は普通に計算できるが、右辺第一項を計算するにはその時刻における各粒子位置での圧力 $P _i$ が必要である。そこで今回は半隂解法を用いる。すなわち、はじめに粘性項・重力項を用いて粒子速度を更新し、その情報を用いて圧力を計算、その後圧力項による粒子速度の計算を行うことにする。すなわち次のよう。

$$
\begin{aligned}
& \dfrac{\boldsymbol{u}^{k+\frac{1}{2}} _i-\boldsymbol{u}^k _i}{\Delta t} = \nu \left(\nabla^2 \boldsymbol{u}\right)^{k} _i + \boldsymbol{g}\newline
& \dfrac{\boldsymbol{u}^{k+1} _i-\boldsymbol{u}^{k+\frac{1}{2}} _i}{\Delta t} = -\dfrac{1}{\rho^0}\left(\nabla P\right)^{k+\frac{1}{2}} _i
\end{aligned}
$$

ただし $k$ は時間方向、 $i$ は空間方向の離散化を表している。その上で、 $\left(\nabla P\right)^{k+\frac{1}{2}} _i$ は次のように導かれるポアソン方程式を用いて求めることができる。まず下式の発散をとると、非圧縮条件より最終的に更新で得られる粒子速度 $\boldsymbol{u}^{k+1} _i$ は $\nabla \cdot \boldsymbol{u}^{k+1} _i = 0$ を満たさなければならないので

$$
\dfrac{-\left(\nabla\cdot \boldsymbol{u}\right)^{k+\frac{1}{2}} _i}{\Delta t} = -\dfrac{1}{\rho^0}\left(\nabla^2 P\right)^{k+\frac{1}{2}} _i
$$

を得る。また、粘性項および重力項を用いて更新した粒子速度 $\boldsymbol{u}^{k+\frac{1}{2}} _i$ は計算途中で現れる量であるが、その計算過程においても連続の式は依然として成立する。すなわち、

$$
\dfrac{D\rho _i}{Dt} + \nabla \cdot\left( \rho _i \boldsymbol{u}^{k+\frac{1}{2}} _i\right)=0
$$

は成立しているはずである。よって $\rho^k _i = \rho^0 = \text{const}$ に注意しながら常識を離散化すれば、

$$
\dfrac{\rho^{k+\frac{1}{2}} _i-\rho^0}{\Delta t} + \rho^0 \left(\nabla \cdot \boldsymbol{u}\right)^{k+\frac{1}{2}} _i = 0 \; \leadsto \left(\nabla \cdot \boldsymbol{u}\right)^{k+\frac{1}{2}} _i = -\dfrac{1}{\Delta t}\dfrac{\rho^{k+\frac{1}{2}} _i-\rho^0}{\rho^0}\sim -\dfrac{1}{\Delta t}\dfrac{n^{k+\frac{1}{2}} _i-n^0}{n^0}
$$

となる。ただし途中で密度 $\rho _i$ を数密度 $n _i$ によって近似した。よって得られた式を組み合わせれば、圧力について次のポアソン方程式が成立する。

$$
-\dfrac{1}{\rho _0}\left(\nabla P^2\right)^{k+\frac{1}{2}} _i=\dfrac{1}{\Delta t^2}\dfrac{n _i^{k+\frac{1}{2}} -n^0}{n^0}
$$

よって先に示したラプラシアンの離散化を用いれば、次のような逆行列演算として圧力を計算できることがわかる。ただし以降は $\phi^{k+\frac{1}{2}}$ を省略して $\phi^*$ と書くことにする。

$$
-\dfrac{1}{\rho _0}\dfrac{2d}{\lambda^0 n^0}\displaystyle\sum _{j\neq i}\left(P _j^{\ast}-P _i^{\ast}\right)w(|\boldsymbol{r} _j^\ast -\boldsymbol{r} _i^\ast |)=\dfrac{1}{\Delta t^2}\dfrac{n _i^\ast -n^0}{n^0}
$$

### 衝突判定

上では割愛していたが、重力項と粘性項により粒子を移動させたら、圧力項を計算する前に一旦粒子の衝突判定を行い衝突条件を満たすような粒子を反発させる。すなわち、二つの粒子間距離に関する条件 $\left|\boldsymbol{X} _j-\boldsymbol{X} _i\right|\lt　r _\text{col}$ をみたし、かつ粒子が衝突する向きに動いている、すなわち $u _{ij}=\left(\boldsymbol{v} _j-\boldsymbol{v} _i\right)\cdot \boldsymbol{e} _{ij}\lt　0$ が満たされる時、衝突インパクト（力積） $J _{ij}(=F\Delta)$ を反発係数 $e$ のもとに計算し、速度を $\boldsymbol{v} _i' = \boldsymbol{v} _i-\sum _{j}J _{ij}\boldsymbol{n}^\perp _{ij}/m _i$ によって更新する。（その後粒子の位置も更新する。）

具体的な式としては、初め

$$
\boldsymbol{e} _{ij}=\dfrac{\boldsymbol{X} _j-\boldsymbol{X} _i}{\left|\boldsymbol{X} _j-\boldsymbol{X} _i\right|}
$$

として計算し、

$$
\begin{aligned}
&\text{if}\; \left|\boldsymbol{X} _j-\boldsymbol{X} _i\right|\lt r _\text{col}\; \text{and}\; u _{ij}=\left(\boldsymbol{v} _j-\boldsymbol{v} _i\right)\cdot \boldsymbol{e} _{ij}\lt 0\;;\quad && J _{ij} = \left(1+e\right)\dfrac{m _im _j}{m _i+m _j}(-u _{ij})\newline
& \text{otherwise}\quad  && J _{ij} = 0 \newline
\end{aligned}
$$

として衝突インパルスを求め、


$$
\boldsymbol{v}' _{i} =\boldsymbol{v} _{i}-\displaystyle\sum _{j}\dfrac{J _{ij}}{m _i}\boldsymbol{e} _{ij}\quad ,\quad \boldsymbol{X}' _i = \boldsymbol{X} _i + \left(\boldsymbol{v}' _{i} -\boldsymbol{v} _{i}\right)\Delta t
$$

として位置を修正する。


### 圧力項の計算
圧力 $P$ はポアソン方程式によって計算され、その離散化は次のようになるのであった。

$$
-\dfrac{1}{\rho _0}\dfrac{2d}{\lambda^0 n^0}\displaystyle\sum _{j\neq i}\left(P _j^{\ast}-P _i^{\ast}\right)w(|\boldsymbol{r} _j^\ast -\boldsymbol{r} _i^\ast |)=\dfrac{1}{\Delta t^2}\dfrac{n _i^\ast -n^0}{n^0}
$$

この時、粒子 $\boldsymbol{X} _i$ 近傍の粒子数密度 $n^* _i$ をあらかじめ計算する必要がある。また、圧力のポアソン方程式を特にあたり境界条件を定める必要があるが、これは計算された粒子数密度 $n^* _i$ が一定値を下回った時にその粒子を自由表面粒子と判定し、その点において $p _i=0$ というディリクレ境界条件を入れることによって問題が解かれる。

具体的には、

$$
n^* _i = \displaystyle\sum _{j\neq i}w\left(\left|\boldsymbol{X} _i-\boldsymbol{X} _j\right|\right)
$$

によって粒子数密度を計算し、

$$
n _i^* \leq\beta  n^0
$$

の時粒子を自由表面粒子と判定することにする。また、壁面における境界条件については圧力勾配が $0$ となるノイマン境界条件を用いるが、これは壁の外側に圧力値を持たないダミー粒子を配置することで近似的に表現できる。

#### ポアソン方程式

上で与えられたポアソン方程式は次のように書き換えることができる。

$$
\displaystyle\sum _{j\neq i}c _{ij}\left(P _i^{\ast}-P _j^{\ast}\right)=b _i
$$

$$
\left(\text{ただし}\; c _{ij} = \dfrac{w(|\boldsymbol{r} _i^* -\boldsymbol{r} _j^* |)}{\rho _0}\dfrac{2d}{\lambda^0 n^0}\;,\; b _i = \dfrac{1}{\Delta t^2}\dfrac{n _i^* -n^0}{n^0}\right)
$$

よって $P _i = P _i^\ast$ としてベクトル $\boldsymbol{P}$ を定め、行列 $A$ を

$$
A _{ij} = \left\lbrace\begin{aligned}
&\; -c _{ij} \quad && \left(j\neq i\right)\newline
&\; \displaystyle\sum _{j'\neq i} c _{ij'} &&\left(j= i\right)
\end{aligned}\right.
$$

によって定めれば、ポアソン方程式は行列方程式 $A\boldsymbol{P}=\boldsymbol{b}$ によって計算することができ、特に $A$ が対称行列となることを利用すれば共役勾配法などで解くことができる。ただし、実際の計算の上では安定性の向上及び収束性の向上のため、 $A\boldsymbol{P}=\boldsymbol{b}$ を次のように修正する。

$$
\left(A+\dfrac{\kappa}{\Delta t^2}I\right)\boldsymbol{P} = \gamma \boldsymbol{b}
$$

ただし上記の修正及び $\kappa,\beta$ の具体的な選び方は「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」を参考にしている。また、ディリクレ境界条件については、表面粒子として判定された $i$ に対して $P _i = b _i=0$ として計算すれば良い。

#### 共役勾配法を用いた行列ソルバ

上の式で出てくる行列は実対称正定値となっている。対称性はいいとして、正定値性については

$$
\boldsymbol{P}^TA\boldsymbol{P} = \cdots = \dfrac{1}{2}\displaystyle\sum _{i,j}c _{ij}\left(P _i-P _j\right)^2 +  \dfrac{\kappa}{\Delta t^2}P _i^2
$$

から従う。これを説明すると次のよう。

今式の形から半正定値性 $\boldsymbol{P}^TA\boldsymbol{P} \geq 0$ は明らか。そこで $\boldsymbol{P} ^TA\boldsymbol{P} =0$ となる条件について考えるが、まず右辺第一項については $c _{ij}\propto w(|\boldsymbol{r} _i^* -\boldsymbol{r} _j^* |)$ であったのでこれが $0$ となるには「（重み関数が非ゼロとなるような）近傍粒子で圧力値が等しくなる」ことが必要となる。よって近傍粒子を辿っていけば自由表面粒子にたどり着くことを考えると、右辺第一項が $0$ となるためには任意の粒子において $P _i=0$ とならなければならないことがわかる。（ただし、「近傍粒子を辿っていっても自由表面粒子にたりつかないような内部粒子」、すなわち孤立して存在するが数密度があまり小さくない粒子（孤立した二つの粒子が非常に接近してる場合など）が存在する場合においてはこれら粒子に対し $P _i=\text{const}\neq0$ となっても右辺第一項が $0$ となってしまう。ただしこの場合においては右辺第二項が非ゼロとなり、 $\boldsymbol{P}^TA\boldsymbol{P} > 0$ となる。）よって $\boldsymbol{P}^TA\boldsymbol{P} \geq 0\Rightarrow \boldsymbol{P}=0$ なので、確かに $A$ は正定値である。

（※ $\kappa/\Delta t$ の項が正定値性を助けているが、それでも先に述べたように「近傍粒子を辿っていっても自由表面粒子にたりつかないような内部粒子」が存在してしまうと正定値性は弱くなってしまう。そこで 「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」においては、このような粒子に対し対角項を二倍にするなどの処理を行うことで計算を安定化させている。（余裕があればこの処理を入れたい...））

また、 $A$ の係数は近傍粒子以外では $0$ となり疎行列となる。よって理想的には近傍粒子リストを作成し、それらに対してのみ行列の係数の計算を行うのが良い。（コード全体に言えることではある。）ただ、今回はそこまでの余裕はないので普通に行列を作って計算することにする。といっても、行列をそのまま作ると粒子数 $^2$ のメモリを保持することになり現実的ではないので、今回は共役勾配法の各ステップごとに行列の係数を計算して（メモリには保持せず）ベクトルとの積を計算することにする。

また、今 $A\boldsymbol{p}=\boldsymbol{b}$ という方程式に対し、共役勾配法のアルゴリズムは次のよう。

$$
\begin{aligned}
&\boldsymbol{p} _0=0\newline
&\boldsymbol{r} _0 = \boldsymbol{b}-A\boldsymbol{P} _0\newline
&\boldsymbol{d} _0 = \boldsymbol{r} _0\newline
&\text{do}\; l= 0,1,2,\ldots\newline
&\qquad \alpha _l = \dfrac{\left<\boldsymbol{r} _l,\boldsymbol{r} _l\right>}{\left<\boldsymbol{d} _l,A\boldsymbol{d} _l\right>}\newline
&\qquad \boldsymbol{P} _{l+1} = \boldsymbol{P} _l + \alpha _l\boldsymbol{d} _l\newline
&\qquad \boldsymbol{r} _{l+1} = \boldsymbol{r} _{l}-\alpha _l A\boldsymbol{d} _l\newline
&\qquad \text{if}\; \left(\dfrac{\left|\boldsymbol{r} _l\right|}{\left|\boldsymbol{b}\right|}<\varepsilon\right)\; \text{break}\newline
&\qquad \beta _l = \dfrac{\left<\boldsymbol{r} _{l+1},\boldsymbol{r} _{l+1}\right>}{\left<\boldsymbol{r} _l,\boldsymbol{r} _l\right>}\newline
&\qquad \boldsymbol{d} _{l+1} = \boldsymbol{r} _{l+1}+\beta _l\boldsymbol{d} _l\newline
&\text{end do}
\end{aligned}
$$

## 陽解法

半陰解法の難点は、計算量の多さである。今回は近傍粒子探査を実装していないので基本的に $\mathcal{O}(N^2)$ の計算量となているが、実装上は近傍粒子探索を入れることにより基本的な計算量が $\mathcal{O}(N)$ 、共役勾配法を用いた圧力計算部分だけ $\mathcal{O}(1.5N)$ となり、実質的にはポアソン方程式の球解が計算のボトルネックとなる。そこで、陰的に圧力を計算するのを辞め、陽的に圧力を求めることを考える。具体的には、半陰解法と同様重力項と粘性項のみを用いて粒子位置を更新したあと、

$$
P^{*} _i = c^2\dfrac{\rho^0}{n^0}\left(n _i^*-n^0\right) 
$$

によって圧力を計算するという考え方である。ただし $c$ は流体の音速である。上式は、音速の定義 $c^2=\partial  P/\partial \rho $ に従って圧力を $P^* _i\sim P^0+c^2 \left(\rho^* _i - \rho^0\right)$ テイラー展開し、 $\rho _i^* \sim \rho _0\frac{n _i^*}{n^0}$ 、 $P^0=0$ と考えることによって得られる。

また、勾配計算についても半陰解法とは異なった方法を用いることになる。今半陰解法においては

$$
\left(\nabla P\right)^* _i = \dfrac{d}{n^0}\displaystyle\sum _{i\neq j}\dfrac{P^* _j-P^* _i}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|^2}\left(\boldsymbol{X} _{j}-\boldsymbol{X} _{i}\right)w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

のように圧力勾配を計算するが、陽解法においては

$$
\left(\nabla P\right)^* _i=\dfrac{d}{n^0}\displaystyle\sum _{i\neq j}\dfrac{P^* _j+P^* _i}{|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|^2}\left(\boldsymbol{X} _{j}-\boldsymbol{X} _{i}\right)w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

のように差分を和で置き換えて計算することになる。詳しくは「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」を参照。また重み関数も陽解法とは異なったものを用いることが望ましく、具体的には

$$
w\left(r\right) = \dfrac{r _e}{r}+\dfrac{r}{r _e}-2
$$

とする。ただし圧力勾配の計算においては

$$
w _\text{grad}\left(r\right) = \dfrac{r _e}{r}-\dfrac{r}{r _e}
$$

を用いることにする。

### 音速の決定

陽解法を用いるにあたっては、流体音速 $c^2$ を定めなkればならない。本来は物理的に正しい値を用いるのが理想的だが、その場合クーラン条件（安定性条件）を満たすために $\Delta t$ を相当小さくしなければならず、計算コストの点から好ましくない。具体的には、今数値安定性条件は上記の音速で記述されるものと、龍朔kの上限値で記述されるものの二つが存在していて次のようになっている。

$$
\dfrac{u _\text{max}\Delta t}{l^0} < 0.2
$$

$$
\dfrac{c\Delta t}{l^0} < 1.0
$$

半陰解法においては常しきだけを満たせば数値安定性を達成できたため $\Delta t$ を相対的にそこまで細かく切る必要はなかったが、陽解法では下式も同時に満たす必要があり、これがネックとなって $\Delta t$ に厳しい制約がかかってしまう。これを回避する一つの方法として、音速を物理的な値ではなく便宜的に $c = u _\text{max}/0.2$ と置く方法が「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」では紹介されている。

この場合、上の数値安定性条件を満たすような半陰解法と同様の時間ステップ幅で親しきの数値安定性条件も満たされるため、計算コストの問題がなくなる。ただしこれは音速を本来の物理的な値よりも小さく見積もっていることになり、実質的に圧縮性を大きく見積もっていることに相当する。すなわち非圧縮条件が崩れ、解が非圧縮解からズレることになる。すなわち、陽解法によって生じる計算コストの問題の解消は「回の精度の悪化」という別の問題へと移し替えることによって解消される（？）ことになる。


<br><br><br>


***
# プログラムの概要
***
<br>

実際のプログラムで行われる処理（サブルーチン）を下に示す。

```fortran:
program main

|-water _tank _and _water _column _2d
|-calConstantParameter
|-mainLoopOfSimulation   
  |-calGravity          
  |-calViscosity        
  |-moveParticle        
  |-collision
  |-calPressure
    |-calNumberDensity
    |-setBoundaryCondition
    |-setSourceTerm
    |-calculate _pressure _by _cg
        |-apply _matrix _to _vector
        |-inner _product _of _vector
    |-removeNegativePressure
    |-setMinimumPressure
  |-calPressureGradient
  |-moveParticleUsingPressureGradient
  |-if (timestep = outputstep) {writeData _inVtuFormat}
  |-if (time>finish time){exit mainloopOfSimulation}

end program
```

## それぞれのsubroutineの説明

### water _tank _and _water _column _2d
　初期の粒子配置を設定するルーチン。
```txt
|            |                                                     |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            | _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|           | 
|                                                                              | 
|                                                                              | 
|                                                                              |  
| _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 
```
こんな感じで水槽と水柱を配置する。具体的には、必要になる粒子の数を```numberofparticles```として、i=1~numberofparticlesの各粒子に対して、粒子の位置を表す```particleposition(i)```や粒子の種類(壁など)を表す```particletype(i)```などを、メモリを確保しながら（allocateしながら）入力していく。

### calConstantParameter

$n^0$ や $\lambda^0$ を計算する。

### mainLoopOfSimuation

k〜k+1ステップに粒子の情報を更新するルーチン。以下のcalGravity~movePariceUsingPressureGradientで構成される。

改めてナビエストークス方程式は

$$\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}$$

のようになっており、mainLoopOfSimuationにおいては、まず重力項と粘性項をcalGravity、calViscosityによって計算し、それら加速度を用いて一旦粒子の情報をアップデートする。次にcalPressureを用いて粒子の圧力を計算し、圧力項をcalPressureGradientによって計算する。その後、それを用いて再度粒子の情報をアップデートする。

#### calGravity

重力項による粒子の加速を計算するルーチン

#### calViscosity

粘性項による粒子の加速を計算するルーチン

#### moveParticle

calGravity,calViscosityで計算した加速度をつかって粒子の位置・速度を更新するルーチン

#### collision

例外処理。異常接近した粒子があった場合に、粒子間距離を広げる。

#### calPressure

各粒子位置での圧力を計算するルーチン。以下のcalNumberDensity~setMinimumPressureで構成される。今、圧力$P$はポアソン方程式を解くことによって得られ、その離散化は

$$
-\dfrac{1}{\rho _0}\dfrac{2d}{\lambda^0 n^0}\displaystyle\sum _{j\neq i}\left(P _j^{k+1}-P _i^{k+1}\right)w(|\boldsymbol{r} _j^* -\boldsymbol{r} _i^* |)=\dfrac{1}{\Delta t^2}\dfrac{n _i^* -n^0}{n^0}
$$

と表されるのであった。よって、圧力を計算するためには粒子近傍の粒子数密度 $n^* _i$ を計算し、その値によって自由表面境界にある粒子を選別した上で自由表面における境界条件 $p _i=0$ を与えなければならない。また、その後のポアソン方程式を解く過程も、「右辺（ソース項の決定）」「左辺の圧力行列の準備」「行列方程式の求解による各粒子点での圧力の計算」「負圧の除去」といったプロセスを辿ることになる。また、後々圧力項（ $\nabla p$ ）を計算するにあたって、近傍最小圧力を計算しておく必要がある。

よって calPressureの中身は次のようになる。

- <ins>calnumberDensity</ins>：粒子数密度を計算するルーチン
- <ins>setBoundaryCondition</ins>：自由表面境界を探すルーチン
- <ins>setSourceTerm</ins>：ポアソン方程式の右辺を準備するルーチン
- <ins>setMatrix</ins>：ポアソン方程式の左辺を準備するルーチン
- <ins>calculate _pressure _by _cg</ins>：ガウスの消去法によるポアソン方程式の求解
  - <ins>apply _matrix _to _vector</ins>：行列の係数を都度計算しベクトルに作用させるルーチン
  - <ins>inner _product _of _vector</ins>：ガウスの消去法で表れる内積計算を、特に近傍粒子に限定して行うルーチン
- <ins>setNegativePressure</ins>：負圧の除去
- <ins>setMinimumPressure</ins>：近傍最小圧力の計算

#### calPressureGradient

圧力項による粒子の加速を計算するルーチン

#### moveParticleusingPressureGradient

calPressureGradientで計算した加速度を使って粒子の位置・速度を更新するルーチン。ついでに、加速度をゼロにリセットする。


<br><br><br>

***
# 各ルーチンの詳細について
***

<br>

## water _tank _and _column _2d
２次元の水柱に関しての初期状態を定めるルーチン。ルーチンの概要は以下のよう。

```fortran
subroutine water _tank _and _water _column _2d(x _watertank,y _watertank,x _watercolumn,y _watercolumn,wallthickness,dummywallthickness)

    nx _watertank = ...
    number _of _particles = ...
    allocate(particle _position(number _of _particles,3))
        .
        .
        .
    
    do iY = ...
        do iX =...
        !particle _position(i,:)=...
        end do 
    end do

end subroutine
```

まず、引数として
```fortran
x _watertank,y _watertank,x _watercolumn,y _watercolumn,wallthickness,dummywallthickness
```
の6つを受け入れる。```x _watertank,y _watertank```で水槽の（内壁の）大きさを定め、```x _watercolumn,y _watercolumn```で水柱の大きさを定め、```wallthickness,dummywallthickness```で各種壁の厚みを定める。いずれも受け入れ引数の単位は[m]。
　で、あらかじめ設定しておいた初期粒子間距離[m]を用いて、トータルで必要な粒子数```numberofparticles```を求め、```particle _position```等に必要粒子数分のメモリを確保していく。
 　その後、```iX,iY```によるループを用いて

$$
(\text{iX}\times\text{初期粒子間距離},\text{iY}\times\text{初期粒子間距離},)
$$

に位置する粒子についての```particletype(i)```などの情報を入れていく。このループは並列化して行う。

## calConstParameter
計算で用いる定数を計算するサブルーチン。まず、

$$
n^0 = \displaystyle\sum _{i\neq j}w(|\boldsymbol{X} _j-\boldsymbol{X} _i|)
$$

の値を計算する。このとき、用いる影響半径 $R _e$ は $n^0$ の用途によって異なることに注意する。

例えばラプラシアン計算用の $n^0$ とついでに $\lambda _0$ を計算するプログラムの概要は、
```fortran
subroutine calc _n0 _for _lambda _and _lambda _for _laplacian


    i _for _Re = floor(Re _for _laplacian/particle _distance)+1
        
    if (dimention == 2) then
    
        do iX = -i _for _Re,i _for _Re  
            do iY = -i _for _Re,i _for _Re  

                xdist = real(iX,8)*particle _distance
                ydist = …

                distance = sqrt(xdist**2 + ydist**2)

                if(iX == 0 .nor. iY == 0) then
                    w = weight _function(distance,Re _for _laplacian) 
                end if

                n0 _for _laplacian = n0 _for _laplacian + w
                lambda = lambda + (distance**2) * w
                
            end do
        end do
        lambda _0 = lambda/n0 _for _laplacian
    end if


end subroutine calc _n0 _and _lambda _for _laplacian

```
みたいな形。

```fortran
i _for _Re = floor(Re _for _laplacian/particle _distance)+1
```

で、初期粒子数密度 $n^0$ の計算に用いる液体内部粒子のおおよその範囲を定めている。この範囲内で、初期粒子間距離に従って粒子を配置して、諸々の値を計算しているイメージ。



## calGravity

重力による加速を求めるルーチン。

$$
\left.\dfrac{D\boldsymbol{u}}{Dt}\right| _{gravity} = g\cdot\hat{\boldsymbol{e}} _y
$$

であるから、
```fortran
Gy=-9.80665

acceleration(i,1) = 0
acceleration(i,2) = Gy
acceleration(i,3) = 0
```
のように計算すれば良い。

## calViscosity

改めて粘性項の離散化は、先に示したラプラシアンの離散化を用いて

$$
\nabla^2 \boldsymbol{u}= \dfrac{\nu \cdot 2d}{\lambda _0 n^0}\displaystyle\sum _{i\neq j}(\boldsymbol{u}(\boldsymbol{X} _{j})-\boldsymbol{u}(\boldsymbol{X} _{i}))w(|\boldsymbol{X} _{j}-\boldsymbol{X} _{i}|)
$$

となる。動粘性係数（単位：m^2/s）については、適当な文献を参照すればいい。これは自由に変更できるパラメタになっている。概要は以下。

```fortran
subroutine calviscosity()

    a = (viscosity*2.0 _8*dimention)/(n0 _for _laplacian*lambda _0)

    loopi : do i = 1,number _of _particles

        loopj : do j=1,number _of _particles

            xdist = particle _position(j,1)-particle _position(i,1)
            ydist = ...

            distance = sqrt(xdist**2.0 _8 + ydist**2.0 _8 + zdist**2.0 _8)

            if (distance < Re _for _laplacian*1.1 _8) then
                w = weight _function(distance,Re _for _laplacian)
                viscosity _term(1) = viscosity _term(1)+(particle _velocity(j,1)-particle _velocity(i,1))*w
                viscosity _term(2) = ...
            end if

        end do loopj

        viscosity _term(:) = a*viscosity _term(:)

        particle _acceleration(i,1)=particle _acceleration(i,1)+viscosity _term(1)
        particle _acceleration(i,2)=...

    end do loopi

end subroutine calviscosity
```

## moveparticle

計算された加速度に従って粒子の位置を更新する。概要は以下。

```fortran
subroutine moveparticle()

    do i= 1,number _of _particles
        if (particle _type(i) == fluid) then

            particle _velocity(i,1) = particle _velocity(i,1)+particle _acceleration(i,1)*time _interval
            particle _velocity(i,2) = ...

            particle _position(i,1) = particle _position(i,1)+particle _velocity(i,1)*time _interval
            particle _position(i,2) = ...
        end if
    end do

end subroutine moveparticle

```

## collision

計算は先に述べた通り以下のように行う。


$$
\boldsymbol{e} _{ij}=\dfrac{\boldsymbol{X} _j-\boldsymbol{X} _i}{\left|\boldsymbol{X} _j-\boldsymbol{X} _i\right|}
$$

$$
\begin{aligned}
\text{if}\; \left|\boldsymbol{X} _j-\boldsymbol{X} _i\right|<r _\text{col}\; \text{and}\; u _{ij}=\left(\boldsymbol{v} _j-\boldsymbol{v} _i\right)\cdot \boldsymbol{e} _{ij}<0\;;\quad & J _{ij} = \left(1+e\right)\dfrac{m _im _j}{m _i+m _j}(-u _{ij})\newline
 \text{otherwise}\; ; \; \quad  & J _{ij} = 0 \newline
\end{aligned}
$$

$$
\boldsymbol{v}' _{i} =\boldsymbol{v} _{i}-\displaystyle\sum _{j}\dfrac{J _{ij}}{m _i}\boldsymbol{e} _{ij}\quad ,\quad \boldsymbol{X}' _i = \boldsymbol{X} _i + \left(\boldsymbol{v}' _{i} -\boldsymbol{v} _{i}\right)\Delta t
$$

ただし、 $m _i = m _j$ と考えるので、結局 $\frac{m _im _j}{m _i+m _j}\frac{1}{m _i} = \frac{1}{2}$ となる。なので $m _i,m _j$ は計算する必要がない。概要は以下。

```fortran
subroutine collision

    do i= 1,number _of _particles
    if (particle _type(i) == fluid) then
        velocity _ix = particle _velocity(i,1)
        velocity _iy = ...
        velocity _after _collision(i,1) = velocity _ix
        velocity _after _collision(i,2) = ...

        do j=1,number _of _particles

            xij = particle _position(j,1)- particle _position(i,1)
            yij = ...
            distance2 = (xij*xij + yij*yij + zij*zij) 

            if (distance2 < collision _distance2 .and. distance2 > 0.0) then
                distance = sqrt(distance2)
                relative _speed _negative = (velocity _ix-particle _velocity(j,1))*xij/distance  + ...

                if (relative _speed _negative > 0.0) then 
                    forceDT = ((1.0+restitution _coefficient)/2.0)*relative _speed _negative
                    velocity _ix = velocity _ix - forceDT*xij/distance
                    velocity _iy = ...

                end if
            end if
        end do

        velocity _after _collision(i,1) = velocity _ix
        velocity _after _collision(i,2) = ...

    end if
    end do

    do i= 1,number _of _particles
        if (particle _type(i) == fluid) then
            particle _position(i,1) = particle _position(i,1) + (velocity _after _collision(i,1)-particle _velocity(i,1))*time _interval
            particle _position(i,2) = ...
            particle _velocity(i,1) = velocity _after _collision(i,1)
            particle _velocity(i,2) = ...
        end if
    end do

end subroutine collision
```

## calnumberDensity

粒子近傍の粒子数密度を計算するルーチン。ただし、粒子数密度は重み関数を用いて

$$
n^* _i = \displaystyle\sum _{j\neq i}w\left(\left|\boldsymbol{X} _i-\boldsymbol{X} _j\right|\right)
$$

のように計算する。

## setboundarycondition

粒子が自由表面粒子であるか否かを判定するルーチン。粒子数密度が一定値を下回った時に、その粒子は流体表面にあると判定する。具体的には

$$
n _i^* \leq\beta  n^0
$$

を満たす粒子を自由表面粒子と判定する。コードでは ```threshhold _ratio _of _number _density``` が $\beta$ に対応している。

## setsourceterm

## writeData _inVtuFormat
VTK（.vtu）フォーマットに出力するためのサブルーチン。VTKフォーマットは可視化のためのフォーマット。空行と空白を同様に扱う。VTKファイルは基本は以下の様な構造になっている.

```xml
<?xml version='1.0',encoding='UTF-8'?>
<VTKFile xmlns='VTK' byte _order='LittleEndian' version='0.1' type='UnstructuredGrid'>
    <unstructuredGrid>
        <Piece NUmberOfcells='50000' NumberOfPoints='50000'>

            <Points>
                <Dataarray NumberOfComponents='2or3' type='Float32' Name='Position' format='ascii'>
                    <!座標データ、n行2/3列、nは粒子数、2/3は空間の次元>
                </Dataarray>
            </Points>

            <PointData>
                <Dataarray NumberOfComponents='i' type='Int32/Float32/etc' Name='hoge' format='ascii'>
                    <!粒子のhoge属性に関するデータ、n行i列、nは粒子数iはデータの次元>
                </Dataarray>
                <Dataarray NumberOfComponents='j' type='Int32/Float32/etc' Name='Foo' format='ascii'>
                    <!粒子のFoo属性に関するデータ、n行j列、nは粒子数iはデータの次元>
                </Dataarray>
                    .
                    .
                    .
            </PointData>

            <Cells>
                <Dataarray　type='Int32' Name='connectivity' format='ascii'>
                    <!後ほど解説>
                </Dataarray>
                <Dataarray　type='Int32' Name='offsets' format='ascii'>
                    <!後ほど解説>
                </Dataarray>
                <Dataarray type='Int32' Name='types 'format='ascii'>
                    <!後ほど解説>
                </Dataarray>
            </Cells>
            
            <Celldata>
                <Dataarray NumberOfComponents='i' type='Int32/Float32/etc' Name='hage' format='ascii'>
                    <!cellのhage属性に関するデータ、m行i列、mはcellの数、iはデータの次元>
                </Dataarray>
                <Dataarray NumberOfComponents='j' type='Int32/Float32/etc' Name='Fooo' format='ascii'>
                    <!cellのFooo属性に関するデータ、m行j列、mはcellの数、jはデータの次元>
                </Dataarray>
                    .
                    .
                    .
            </Celldata>
        </Piece>
    </unstructuredGrid>
</VTKFile>
```
主要な構造だけ抜き出すと以下のよう。[参考](https://docs.vtk.org/en/latest/vtk _file _formats/vtkxml _file _format.html?utm _source=chatgpt.com)
```xml
<VTKFile type="UnstructuredGrid" ... >
  <UnstructuredGrid>
    <Piece NumberOfPoints="#" NumberOfCells="#">
        <Points>...<!-座標->...</Points>
        <PointData>...<!-あ->...</PointData>
        <Cells>...</Cells>
        <CellData>...</CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>

```
<br>

### VTKFile 条件子
```xml
<VTKFile xmlns='VTK' byte _order='LittleEndian' version='1.0' type='UnstructuredGrid'>
```
で指定される、VTKファイルに関する諸々のパラメタ。

<ins>byte _order</ins>：エンディアン。（データを前から並べるか、後ろから並べるか、くらいの認識。）macはLittleEndian。

<ins>version</ins>：VTKフォーマットの仕様のバージョン。1.0でいいのかな。

<ins>type</ins>：stucturedは格子データのような規則正しく並んだデータ、unstructuredは粒子データのような不規則に並んだデータ。色々あるっぽいが、粒子法では```unstructuredGrid```を用いる。unstructuredGridは拡張子```.vtu```に対応しているので、このコードでは```.vtu```ファイルに出力している。
<br>

### piece条件子
```NumberOfPoints```は粒子数、```NumberOfCell```はセルの数。先に下のpointとcellに関する解説を読んだ方がわかりやすいと思う。
<br>

### points
```xml
<Points>
    <Dataarray NumberOfComponents='2or3' type='Float32' Name='Position' format='ascii'>
        0.0 0.0 0.0
        0.0 0.0 0.1
        0.0 0.0 0.2
            .
            .
            .
    </Dataarray>
</Points>
```
各粒子の座標情報が入る。```NumberOfComponents```は空間の次元。xmlは改行と空白を区別しないので、
```xml
0.0 0.0 0.0
0.0 0.0 0.1
0.0 0.0 0.2
    .
    .
    .
```
の部分は
```xml
0.0 0.0 0.0 0.0 0.0 0.1 ...
```
のように認識される。よって、```NumberOfComponents='3'```で、「この座標は三次元やで。やから３つごとにデータを読み取りなさいな」と言う指定を入れてやる必要がある。```type```はデータの型（INTなど）、```name```はまあ、自由に設定していいそのDataarrayの名前。```Format```はasciiかbinary。大規模計算とかだとbinaryの方がいいらしい。
<br>

### pointdata
各粒子にデータを載せる。圧力だったり温度だったり粒子の種類だったり。例えば温度という一次元データを各粒子に持たせることを考えると、
```xml
<Dataarray NumberOfComponents='1' type='Int32/Float32/etc' Name='temperature' format='ascii'>
    273
    280
    250
     .
     .
     .
</Dataarray>
```
ここでも、```NumberOfComponetnsは次元。paraviewなら三次元データ（速度場とか）も可視化できるそうな。すごいッピ！
### cells
ざっくばらんにいって仕舞えば、メッシュ。一つのメッシュのことをセルと言う。例えば下の i、j、k、l 番目の粒子は四角形のセルを作っている。
```
   i . _ _. k
     |  |
   j . _ _. l
```
そんな感じに、いろんな形のセルをいくらでも作ることができて、それを指定するのが```<cells>```の項目。（つまり、粒子集合の部分集合がセル）必ず、以下の三つをもつ。（それぞれ```Name```固定。）
```xml
<Cells>
    <Dataarray type='Int32' Name='connectivity' format='ascii'>
    <Dataarray type='Int32' Name='offsets' format='ascii'>
    <Dataarray type='Int32' Name='types 'format='ascii'>
</Cells>
```
<br>

<ins>**connectivity**</ins>：まずconnectivityだが、これはcellの頂点に関するデータ。
```xml
<Dataarray type='Int32' Name='connectivity' format='ascii'>
    0 3 5
    3 5 7 8
      .
      .
      .
</Datarray>
```
上のコードなら、1,4,6番目の粒子で三角形のセルを作り、4,6,8,9番目で４角形のセルを作り...という風。（vtkファイルの配列は0始まりなので、```<points>```で定めた粒子において、１番目の粒子のindexが0となる。~~0から始まるの直感的じゃなさすぎてほんと嫌い。Fortranを見習ってほしい。~~）ただ、何度も述べるように、xmlは改行を認識しないので実際には
```xml
<Dataarray type='Int32' Name='connectivity' format='ascii'>
    0 3 5 3 5 7 8...
</Datarray>
```
のように認識されている。よって、どこでcellを区切るのかを明示してやる必要がある。それをするのが以下の```offset```
<br>

<ins>**offset**</ins>
```xml
<Dataarray type='Int32' Name='offsets' format='ascii'>
    3
    7
    .
    .
    .
</Dataarray>
```
上のコードは、```connectivity```において「三番目までの粒子番号が最初のセル、（四番目から）七番目までの粒子番号が二つ目のセル、...」という指定を行う。
<br>

<ins>**types**</ins>：それぞれのcellが三角形なのか四角形なのか、みたいなことを明示してやる必要がある。それぞれの形状とtype番号との対応は[リンク](https://vtk.org/doc/nightly/html/vtkCellType _8h _source.html)にある。今、
```cpp
  VTK _TRIANGLE = 5,
  VTK _QUAD = 9,
```
なので、```types```は以下のように書けば良いとわかる。
```xml
<Dataarray type='UInt8' Name='types 'format='ascii'>
    5
    9
    .
    .
    .
</Dataarray>
```
ちなみに、```types```の種類的に、```UINT8``` (符号なし8ビット)で十分だったりする。

### celldata
pointdataと同様、cellに温度だったり圧力だったりの情報を載せる。
```xml
<Dataarray NumberOfComponents='1' type='Int32/Float32/etc' Name='pressure' format='ascii'>
    10
    15
    .
    .
    .
</Dataarray>
```
cellの数だけデータ数（行数）が存在することになる。
<br>

### 粒子法におけるcellの指定
粒子法においては、cellは（多分）使う必要がない。なので、各点が一つのcellを作るとしてcellデータを記入すれば良い。この時、typesは```VTK _VERTEX=1```なので、1と記入すれば良い。
```hml
<cells>
    <Dataarray type='Int32' Name='connectivity' format='ascii'>
        0
        1
        .
        .
        .
        n-1
    </Dataarray>
    <Dataarray type='Int32' Name='offsets' format='asciss'>
        1
        2
        .
        .
        .
        n
    </Dataarray>
    <Dataarray type='UInt8' Name='types' format='ascii'>
        1
        1
        1
        .
        .
        .
    </Dataarray>        
</cells>
```
