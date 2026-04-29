---
title: "第1章：確率空間・測度論入門"
---

# 第Ⅰ部 第1章：確率・確率変数・分布の基礎

> **ステータス**：🔄 執筆中
> **執筆開始日**：2026-04-28
> **最終更新日**：2026-04-29

---

## この章で答える問い

「確率とは何か」を厳密に定義し、そこから出発して「確率変数」「分布」「期待値」「大数の法則」「中心極限定理」まで、一本の論理の糸でつなぐ。これらはすべて計量経済学の数学的基盤であり、この章を理解せずに先へ進むことはできない。

---

## 0. 数学の準備

本章の内容を理解するために必要な数学——積分・テイラー展開・行列——を最初にまとめて復習する。知っている人は読み飛ばしてよいが、後で「あの公式はどこから来たのか？」と疑問に思ったときに戻ってほしい。

---

### 0.1 リーマン積分（高校数Ⅲの積分）

高校数学で学ぶ積分は**リーマン積分（Riemann integral）**と呼ばれる。その考え方は非常にシンプルだ。

**考え方：$x$ 軸（横軸）を細かく刻む**

関数 $f(x)$ のグラフと $x$ 軸の間の「面積」を求めたい。区間 $[a, b]$ を $n$ 個の小区間 $[x_0, x_1], [x_1, x_2], \ldots, [x_{n-1}, x_n]$ に分割し、各小区間に「高さ $f(x_k^*)$、幅 $\Delta x_k = x_k - x_{k-1}$」の長方形を立てる。

$$\int_a^b f(x)\,dx = \lim_{n \to \infty} \sum_{k=1}^{n} f(x_k^*)\,\Delta x_k$$

分割を無限に細かくしたとき（$\Delta x_k \to 0$）の長方形の面積の合計が、リーマン積分の定義である。

```
【図：リーマン積分のイメージ】

  y
  │    ／￣＼
  │   ／   ＼
  │  ／  ■■■＼
  │ ／■■■│■■│＼
  │■■■│■│■│■■■│
  └──┼──┼──┼──┼── x
     a          b

  各 ■ が「高さ × 幅」の長方形。
  細かく刻むほど正確な面積になる。
```

**計算例（高校の復習）**：$\int_0^2 x^2\,dx = \left[\frac{x^3}{3}\right]_0^2 = \frac{8}{3}$

**よく使う積分公式**：

$$\int x^n\,dx = \frac{x^{n+1}}{n+1} + C \quad (n \neq -1)$$

$$\int e^x\,dx = e^x + C, \quad \int \frac{1}{x}\,dx = \ln|x| + C$$

$$\int_{-\infty}^{\infty} e^{-x^2/2}\,dx = \sqrt{2\pi} \quad \text{（ガウス積分：重要！）}$$

**部分積分の公式**（期待値の計算でよく使う）：

$$\int u\,dv = uv - \int v\,du \quad \Longleftrightarrow \quad \int f(x)g'(x)\,dx = f(x)g(x) - \int f'(x)g(x)\,dx$$

**リーマン積分の限界**：次の関数はリーマン積分できない。

$$f(x) = \begin{cases} 1 & x \text{ が有理数} \\ 0 & x \text{ が無理数} \end{cases}$$

どんなに細かく $x$ 軸を分割しても、各小区間の中には有理数と無理数が混在するため、長方形の「高さ」が定まらない。しかし確率論では、このような「変則的な」関数を積分したい場面が出てくる。

---

### 0.2 ルベーグ積分（測度論の積分）

ルベーグ積分（Lebesgue integral）は、リーマン積分の「$x$ 軸を刻む」という発想を逆転させた積分である。

**考え方：$y$ 軸（縦軸）を細かく刻む**

「高さが $y_k$ 以上になっている $x$ の集合はどのくらい『広い』か」を測り、$y_k \times \text{（その広さ）}$ を足し合わせる。

```
【図：ルベーグ積分のイメージ】

  y
  │    ／￣＼
y₃│・・・・・・・・・・←この高さ以上のxの集合の"幅"を測る
  │   ／   ＼
y₂│・・・・・・・・・・←
  │  ／       ＼
y₁│・・・・・・・・・・・←
  │
  └──────────── x
      ←――→←→←→
      （y₁以上のxの集合）

  y軸を水平にスライスして、
  各スライスの「xの集合の広さ」を足し上げる。
```

$$\int f\,d\mu \approx \sum_{k} y_k \cdot \mu\!\left(\{x : y_k \leq f(x) < y_{k+1}\}\right)$$

ここで $\mu(\cdot)$ は「集合の広さ」を測る**測度（measure）**である。

**ルベーグ積分の威力**：

- リーマン積分できない関数にも積分を定義できる
- 離散確率変数（さいころ、コイン）と連続確率変数（身長、株価）を**同じ式**で扱える
- 「確率変数の期待値」を離散・連続の場合分けなしに書ける：

$$E[X] = \int_\Omega X(\omega)\,dP(\omega)$$

この式の意味は「実験の全結果 $\omega$ について、$X(\omega)$ の値を確率 $P$ の重みで足し合わせる」である。

**2つの積分の関係まとめ**：

| | リーマン積分 | ルベーグ積分 |
|---|---|---|
| 分割の方向 | $x$ 軸（横）を刻む | $y$ 軸（縦）を刻む |
| 使える関数 | 連続関数・区分連続関数 | （測定可能な）あらゆる関数 |
| 使える状況 | 通常の計算 | 確率論・測度論の理論的基盤 |
| 計算結果 | （両方使える場合は）一致 | 一致 |

**本書の方針**：日常的な計算ではリーマン積分（高校数Ⅲ）で十分である。ルベーグ積分は「なぜその公式が成り立つか」の理論的背景として登場する。証明の詳細は付録Cに譲る。

---

### 0.3 ラドン＝ニコディム微分（わかりやすく）

確率論で「連続確率変数の確率密度関数（PDF）がなぜ積分で確率を与えるか」を説明するときに登場する概念が**ラドン＝ニコディム微分（Radon-Nikodym derivative）**である。

難しい言葉だが、考え方は「密度」のアナロジーで理解できる。

**アナロジー：人口密度**

東京の人口が 1,400万人だとする。東京全体の「人口測度」は 1,400万人である。しかし場所によって人口密度は異なる——渋谷は 1km² あたり 2万人かもしれないし、奥多摩は 1km² あたり 10人かもしれない。

$$\text{地域 } A \text{ の人口} = \int_A \underbrace{f(x)}_{\text{人口密度}} \,d\underbrace{\lambda(x)}_{\text{面積}}$$

ここで $f(x)$ が「人口測度をルベーグ測度（面積）で割ったもの」——つまり人口密度——である。これが**ラドン＝ニコディム微分**の直感だ。

**確率への応用**：

確率測度 $P$（「確率を割り当てるルール」）とルベーグ測度 $\lambda$（「通常の長さ・面積の概念」）があるとき、

$$P(A) = \int_A f(x)\,d\lambda(x)$$

を満たす関数 $f$ が存在する場合、$f$ を「$P$ の $\lambda$ に対するラドン＝ニコディム微分」といい、

$$f = \frac{dP}{d\lambda}$$

と書く。これがまさに**確率密度関数（PDF）**の正体である。

> **一言でいえば**：「確率を長さ（面積）で割ったもの」が確率密度関数。人口密度が「人口を面積で割ったもの」であるのと全く同じ構造だ。

---

### 0.4 テイラー展開・マクローリン展開

**テイラー展開（Taylor expansion）**は、「ある点 $a$ の近くで関数 $f(x)$ を多項式で近似する」技法である。中心極限定理の証明で直接使う。

**定義（テイラー展開）**：関数 $f(x)$ が点 $a$ の近くで何度でも微分できるとき、

$$f(x) = f(a) + f'(a)(x-a) + \frac{f''(a)}{2!}(x-a)^2 + \frac{f'''(a)}{3!}(x-a)^3 + \cdots = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!}(x-a)^n$$

ここで $n! = n \times (n-1) \times \cdots \times 2 \times 1$（$n$ の階乗）、$f^{(n)}$ は $n$ 階導関数。

**定義（マクローリン展開）**：$a = 0$ としたテイラー展開を特に**マクローリン展開（Maclaurin expansion）**という。

$$f(x) = f(0) + f'(0)x + \frac{f''(0)}{2!}x^2 + \frac{f'''(0)}{3!}x^3 + \cdots$$

**理解のポイント**：テイラー展開は「$x = a$ の周辺でだけ正確」な近似式だ。$a$ から離れるほど精度は落ちる。しかし $x$ が $a$ に近い範囲では、高次の項（$(x-a)^3, (x-a)^4, \ldots$）は非常に小さいので、2〜3次の項で打ち切っても十分に近似できる。

**重要なマクローリン展開の一覧**：

以下は計量経済学で繰り返し登場する。

$$e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{n=0}^{\infty} \frac{x^n}{n!}$$

$$\ln(1+x) = x - \frac{x^2}{2} + \frac{x^3}{3} - \cdots \quad (|x| < 1)$$

$$(1+x)^n = 1 + nx + \frac{n(n-1)}{2!}x^2 + \cdots \quad (|x| < 1)$$

$$\sin x = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \cdots, \quad \cos x = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \cdots$$

**中心極限定理で使う展開**（後の第14節で詳しく使う）：

$e^x$ のマクローリン展開で $x$ の代わりに小さな量 $t/\sqrt{n}$ を代入すると、

$$e^{t/\sqrt{n}} = 1 + \frac{t}{\sqrt{n}} + \frac{t^2}{2n} + O\!\left(n^{-3/2}\right)$$

さらに $\left(1 + \frac{c}{n}\right)^n \to e^c$ という極限公式（高校数学で $e$ の定義として登場するもの）も使う。

$$\left(1 + \frac{x}{n}\right)^n \to e^x \quad (n \to \infty)$$

**確認**（$c=1$ の場合）：
$$\left(1 + \frac{1}{n}\right)^n \to e \approx 2.718\ldots$$

> **📌 コラム：連続複利と $(1+x/n)^n \to e^x$ の関係**
>
> この極限公式は、金融の**連続複利（continuous compounding）**のアイデアそのものだ。
>
> 元本 1 円を年利 $r$ で、1年を $n$ 回に分けて複利運用したときの1年後の残高は $(1 + r/n)^n$ 。分割回数 $n$ を増やすほど「利息に利息がつく」頻度が上がる。年利 $r = 100\%$（$x = 1$）での具体的な収束を見てみよう。
>
> | 分割回数 $n$ | 計算 | 1年後の残高 |
> |---|---|---|
> | 1（年1回） | $(1 + 1)^1$ | 2.000000 |
> | 2（半年ごと） | $(1 + 0.5)^2$ | 2.250000 |
> | 12（月1回） | $(1 + 1/12)^{12}$ | 2.613035 |
> | 365（毎日） | $(1 + 1/365)^{365}$ | 2.714567 |
> | 8760（毎時間） | $(1 + 1/8760)^{8760}$ | 2.718127 |
> | $\infty$（連続複利） | $\lim (1 + 1/n)^n$ | $e \approx 2.71828\ldots$ |
>
> 「複利を無限に細かく刻んだ極限」が $e$ の定義そのもの。一般に年利 $r$ の連続複利では $e^r$ 倍になる。
>
> **CLT の証明との接続**：第14節の証明ステップ 4 で登場した $\left(1 + \dfrac{t^2/2}{n}\right)^n \to e^{t^2/2}$ は「年利 $t^2/2$ の連続複利」と全く同じ構造だ。$n$ 回分割で少しずつ増える、その極限が指数関数になる——これが CLT の証明でマクローリン展開が現れる理由の本質である。

---

### 0.5 行列の基礎

計量経済学では複数の変数・複数の観測値を同時に扱う。これを効率的に書くために**行列（matrix）**を使う。ここでは必要最低限を復習する。

**行列の定義**：

$m$ 行 $n$ 列の数値の配列を $m \times n$ 行列という。

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

$a_{ij}$ は $i$ 行 $j$ 列の要素。列ベクトルは $n \times 1$ 行列（$n$ 行 1 列）。

**行列の和（同じサイズ同士）**：

$$A + B = \begin{pmatrix} a_{11}+b_{11} & \cdots \\ \vdots & \ddots \end{pmatrix} \quad \text{（要素ごとに足す）}$$

**スカラー倍**：$cA$（全要素を $c$ 倍）

**行列の積**（$A$ が $m \times k$、$B$ が $k \times n$ のとき $AB$ は $m \times n$）：

$$(AB)_{ij} = \sum_{l=1}^{k} a_{il}\,b_{lj} \quad \text{（$A$ の $i$ 行と $B$ の $j$ 列の内積）}$$

> **重要**：$AB \neq BA$ が一般に成立しない（行列の積は順番を変えると値が変わる）。

**転置行列（transpose）**：行と列を入れ替える。$A$ の $(i,j)$ 要素が $A'$（または $A^T$）の $(j,i)$ 要素になる。

$$A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix} \Rightarrow A' = \begin{pmatrix} 1 & 3 & 5 \\ 2 & 4 & 6 \end{pmatrix}$$

性質：$(AB)' = B'A'$（順序が逆転することに注意）

**単位行列（identity matrix）**：対角成分が全て 1、他が全て 0 の正方行列。

$$I_n = \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix}, \quad AI = IA = A$$

**逆行列（inverse matrix）**：$n \times n$ 正方行列 $A$ に対して $AB = BA = I$ を満たす $B$ を $A$ の逆行列といい $A^{-1}$ と書く。逆行列が存在する $\iff$ $A$ が**正則（行列式 $\det(A) \neq 0$）**。

$$2 \times 2 \text{ 行列の逆行列：} \begin{pmatrix} a & b \\ c & d \end{pmatrix}^{-1} = \frac{1}{ad-bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$$

**対称行列（symmetric matrix）**：$A = A'$ を満たす行列。分散共分散行列は必ず対称。

**正定値行列（positive definite matrix）**：全ての $0 \neq \mathbf{v} \in \mathbb{R}^n$ に対して $\mathbf{v}'A\mathbf{v} > 0$ を満たす対称行列。$(X'X)$ の逆行列が存在するための条件として OLS で登場する。

**計量経済学で最重要の行列公式（第7章で証明する）**：

$$\hat{\beta} = (X'X)^{-1}X'y \quad \text{（OLS 推定量）}$$

この公式の意味は今は理解しなくてよい。「こういう行列の式が登場する」という予告として記しておく。

---

## 1. 確率空間の基礎用語

### 1.1 標本空間（Sample Space）

**標本空間**とは「起こりうる全ての結果を集めた集合」のことだ。記号で $\Omega$（大文字のオメガ）と書く。

- コインを1回投げる → $\Omega = \{表, 裏\}$
- さいころを1回振る → $\Omega = \{1, 2, 3, 4, 5, 6\}$
- 1人の人の身長を測る → $\Omega = (0, \infty)$（正の実数全体）
- 株価の日次変化率を測る → $\Omega = (-\infty, \infty) = \mathbb{R}$（実数全体）

「コインを1回投げる」実験なら「表が出た」「裏が出た」の2通りしかない。しかし「身長を測る」実験では、結果は 170cm かもしれないし 170.01cm かもしれない——結果が連続した値をとる。このような場合、標本空間は実数の区間になる。

**実数直線（Real number line）**：$\mathbb{R} = (-\infty, +\infty)$ は全ての実数を並べた直線のこと。負の無限大から正の無限大まで続く数の列。確率論では「連続した値をとる確率変数の標本空間」として登場する。

### 1.2 σ-加法族と事象

$\Omega$ の部分集合（$\Omega$ の一部の要素を集めた集合）を**事象（event）**という。「偶数が出る」は $\{2, 4, 6\} \subseteq \{1,2,3,4,5,6\}$ という事象だ。

しかし $\Omega$ が実数直線のように無限集合の場合、「全ての部分集合に確率を割り当てる」と矛盾が生じることが知られている。そこで「確率を割り当てられる事象の集まり」を **σ-加法族（sigma-algebra）** と呼んで制限する。

σ-加法族 $\mathcal{F}$ の条件：
- $\Omega \in \mathcal{F}$（全体集合は事象）
- $A \in \mathcal{F} \Rightarrow A^c \in \mathcal{F}$（「Aでない」も事象）
- 可算個の事象の和集合も事象

> **直感**：「この実験で議論できる出来事の一覧リスト」が σ-加法族。一覧に載っている出来事の「逆」も「どちらか一方」もリストに入れるべき——という当然の要請。

### 1.3 確率測度とコルモゴロフの公理

**確率測度（probability measure）** $P$ は「事象に0〜1の数（確率）を割り当てる関数」で、次の3公理を満たす。

$$\text{（K1）非負性}：P(A) \geq 0 \qquad \text{（K2）全確率}：P(\Omega) = 1$$

$$\text{（K3）加算加法性}：A_i \cap A_j = \emptyset \Rightarrow P\!\left(\bigcup_{n=1}^{\infty} A_n\right) = \sum_{n=1}^{\infty} P(A_n)$$

3つ組 $(\Omega, \mathcal{F}, P)$ を**確率空間（probability space）**という。

**3公理から導かれる重要な性質**：

- $P(\emptyset) = 0$（何も起きない確率は0）
- $P(A^c) = 1 - P(A)$（余事象）
- $P(A \cup B) = P(A) + P(B) - P(A \cap B)$（加法定理）
- $A \subseteq B \Rightarrow P(A) \leq P(B)$（単調性）

---

## 2. 確率変数（Random Variable）

### 2.1 確率変数とは何か

「確率変数」という言葉は「値がランダムに変わる変数」という印象を与えるが、数学的にはもっと明確に定義できる。

**定義 2.1（確率変数）**：確率空間 $(\Omega, \mathcal{F}, P)$ 上の**確率変数（random variable）**とは、

$$X : \Omega \to \mathbb{R}$$

という関数（標本空間の各要素に実数を対応させる規則）であって、$\{X \leq x\} = \{\omega \in \Omega : X(\omega) \leq x\} \in \mathcal{F}$（確率が割り当てられる）を全ての実数 $x$ に対して満たすものをいう。

> **直感**：確率変数は「実験の結果 $\omega$ を数値 $X(\omega)$ に変換する関数」だ。サイコロの目を「偶数なら1、奇数なら0」に変換する関数も確率変数。

**表記の約束**（重要）：

- 確率変数は**大文字**：$X, Y, Z$
- その実現値（具体的な値）は**小文字**：$x, y, z$
- 「$P(X = x)$」は「確率変数 $X$ が値 $x$ をとる確率」

**例**：コイン投げ。$\Omega = \{表, 裏\}$。

$$X(\omega) = \begin{cases} 1 & \omega = 表 \\ 0 & \omega = 裏 \end{cases}$$

$X$ は確率変数。$P(X = 1) = P(\{表\}) = 1/2$。

---

## 3. 離散確率変数と連続確率変数

確率変数はとりうる値の構造によって2種類に分かれる。

| | 離散確率変数 | 連続確率変数 |
|---|---|---|
| 値の種類 | 有限個または整数のように数えられる | 区間全体（数えられない） |
| 確率の表現 | 確率質量関数（PMF） | 確率密度関数（PDF） |
| 例 | コインの表の回数、さいころの目 | 身長、気温、株価の変化率 |
| 規格化条件 | $\sum_x p(x) = 1$ | $\int_{-\infty}^{\infty} f(x)\,dx = 1$ |

### 3.1 離散確率変数と確率質量関数（PMF）

**PMF = Probability Mass Function = 確率質量関数**

「確率質量関数」とは「離散確率変数の各値の確率を返す関数」である。「質量（mass）」という言葉は「各点に確率の重みが乗っている」イメージからきている。

**定義 3.1（確率質量関数・PMF）**：

$$p(x) = P(X = x), \quad p(x) \geq 0, \quad \sum_{x} p(x) = 1$$

**例（さいころ）**：$p(x) = 1/6$（$x = 1, 2, 3, 4, 5, 6$）、他は 0。

**例（ベルヌーイ分布）**：$P(X=1) = p$、$P(X=0) = 1-p$。これは「成功か失敗か」の実験（コインの表裏、選挙での当選・落選など）を表す。

### 3.2 連続確率変数と確率密度関数（PDF）

**PDF = Probability Density Function = 確率密度関数**

連続確率変数では $P(X = x) = 0$（特定の1点の確率はちょうど0）なので、PMF のような「1点の確率」では表現できない。代わりに「確率の密度（densityは濃さという意味）」という概念を使う。

**定義 3.2（確率密度関数・PDF）**：関数 $f : \mathbb{R} \to [0, \infty)$ が

$$P(a \leq X \leq b) = \int_a^b f(x)\,dx$$

を全ての $a \leq b$ に対して満たすとき、$f$ を $X$ の**確率密度関数（PDF）**という。$f$ は必ず：

$$\text{（D1）}\ f(x) \geq 0 \quad \text{（D2）}\ \int_{-\infty}^{\infty} f(x)\,dx = 1$$

を満たす。

**確率は「面積」である**：

```
  f(x)
   │      ／＼
   │     /  \
   │    / ■■ \   ← この部分の面積 = P(a≤X≤b)
   │   /■■■■ \
   │  /        \
   └───────────── x
         a   b
```

> **重要な注意**：$f(x)$ は「確率」ではなく「確率の濃さ（密度）」である。$f(x) = 3$ でも問題ない（面積が1になれば十分）。確率は必ず積分（面積）で求める。

**ルベーグ積分の視点**：連続確率変数の確率測度 $P$ とルベーグ測度（通常の長さの概念）$\lambda$ の間に

$$P(A) = \int_A f(x)\,d\lambda(x)$$

が成り立つとき、$f = dP/d\lambda$ が「確率密度関数（PDF）」の正体であり、これがラドン＝ニコディム微分（0.3節参照）だ。

---

## 4. 累積分布関数（CDF）

**CDF = Cumulative Distribution Function = 累積分布関数**

「累積（Cumulative）」とは「積み上げ」という意味。「$X$ が $x$ 以下の値をとる確率」を返す関数だ。

**定義 4.1（累積分布関数・CDF）**：

$$F(x) = P(X \leq x) = \begin{cases} \displaystyle\sum_{t \leq x} p(t) & \text{（離散：$x$ 以下の全 $t$ の確率を合計）} \\[8pt] \displaystyle\int_{-\infty}^{x} f(t)\,dt & \text{（連続：$-\infty$ から $x$ まで密度を積分）} \end{cases}$$

**CDF の 4 つの性質**（必ず成立する）：

| 性質 | 内容 | 直感 |
|-----|------|-----|
| 単調非減少 | $x_1 < x_2 \Rightarrow F(x_1) \leq F(x_2)$ | 上限 $x$ が増えるほど確率は増える |
| 右連続 | $\lim_{h \downarrow 0} F(x+h) = F(x)$ | 右から近づく極限が値と一致 |
| 左極限が 0 | $\lim_{x \to -\infty} F(x) = 0$ | $-\infty$ 以下には何もない |
| 右極限が 1 | $\lim_{x \to +\infty} F(x) = 1$ | 何かは必ず起きる |

**よく使う計算式**：

$$P(a < X \leq b) = F(b) - F(a)$$

**PDF と CDF の関係**（連続の場合）：

$$F(x) = \int_{-\infty}^{x} f(t)\,dt \quad \iff \quad f(x) = \frac{d}{dx}F(x)$$

> **なぜ CDF が便利か**：離散・連続どちらの確率変数にも定義でき、「逆関数法」でシミュレーションを生成したり、検定統計量の分布の計算に使ったりする。

---

## 5. 代表的な確率分布

### 5.1 ベルヌーイ分布（Bernoulli Distribution）—— $\text{Ber}(p)$

「成功（1）か失敗（0）か」の2値をとる最も基本的な分布。名前はスイスの数学者ヤーコブ・ベルヌーイ（Jakob Bernoulli）に由来。

$$P(X = 1) = p, \quad P(X = 0) = 1 - p \quad (0 < p < 1)$$

$$E[X] = p, \quad \text{Var}(X) = p(1-p)$$

**例**：コインを1回投げて表が出たか（$p = 0.5$）、薬が効いたか（$p$ は薬の効果）。

### 5.2 二項分布（Binomial Distribution）—— $\text{Bin}(n, p)$

ベルヌーイ試行を $n$ 回繰り返したときの「成功回数」の分布。二項分布の名は二項定理（$(a+b)^n$ の展開）に由来。

$$P(X = k) = \binom{n}{k} p^k (1-p)^{n-k}, \quad k = 0, 1, \ldots, n$$

ここで $\binom{n}{k} = \frac{n!}{k!(n-k)!}$ は二項係数（$n$ 個から $k$ 個を選ぶ組み合わせの数）。

$$E[X] = np, \quad \text{Var}(X) = np(1-p)$$

**例**：コインを 10 回投げたとき表が $k$ 回出る確率（$n=10, p=0.5$）。

### 5.3 一様分布（Uniform Distribution）—— $U(a, b)$

区間 $[a, b]$ 上のどこも等しく起こりやすい分布。「一様（Uniform）」は「均等」を意味する。

$$f(x) = \frac{1}{b-a} \quad (a \leq x \leq b), \quad f(x) = 0 \quad \text{（それ以外）}$$

$$F(x) = \frac{x-a}{b-a} \quad (a \leq x \leq b)$$

$$E[X] = \frac{a+b}{2}, \quad \text{Var}(X) = \frac{(b-a)^2}{12}$$

**例**：乱数生成（コンピュータの乱数は $U(0,1)$ を出発点とする）。

### 5.4 正規分布（Normal Distribution）—— $N(\mu, \sigma^2)$

計量経済学で最も重要な分布。「ガウス分布」とも呼ばれる（カール・フリードリヒ・ガウスに由来）。標準偏差 $\sigma$ を中心に左右対称な「釣り鐘型」をしている。

**標準正規分布** $N(0, 1)$ の PDF：

$$\phi(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}$$

**一般正規分布** $N(\mu, \sigma^2)$ の PDF：

$$f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$$

$$E[X] = \mu, \quad \text{Var}(X) = \sigma^2$$

**CDF（標準正規）**：$\Phi(x) = \int_{-\infty}^x \phi(t)\,dt$（閉じた式なし → 数表・計算機）

$$\Phi(1.96) \approx 0.975, \quad \Phi(1.645) \approx 0.95, \quad \Phi(2.576) \approx 0.995$$

**標準化（standardization）**：$X \sim N(\mu, \sigma^2)$ のとき $Z = \frac{X-\mu}{\sigma} \sim N(0,1)$

**$\int \phi(x)dx = 1$ の証明（ガウス積分）**：

$I = \int_{-\infty}^{\infty} e^{-x^2/2}\,dx$ とおく。$I^2$ を計算する。

$$I^2 = \left(\int_{-\infty}^{\infty} e^{-x^2/2}\,dx\right)\!\left(\int_{-\infty}^{\infty} e^{-y^2/2}\,dy\right) = \int_{-\infty}^{\infty}\!\int_{-\infty}^{\infty} e^{-(x^2+y^2)/2}\,dx\,dy$$

**【図形的イメージ：なぜ極座標変換がうまくいくのか】**

$e^{-(x^2+y^2)/2}$ という関数は、原点からの距離 $r = \sqrt{x^2+y^2}$ だけに依存する。つまり $xy$ 平面上でこの関数の「等高線」を描くと、すべて**原点を中心とする同心円**になる。

```
         y
         ↑
    2 ─  |  ～～～（r=2、高さ e^{-2} ≒ 0.135）
    1 ─  | ～～（r=1、高さ e^{-0.5} ≒ 0.607）
    0 ─  ●──────────────→ x
   -1 ─  | ～～
   -2 ─  |  ～～～
         |
   （● が原点。等高線はすべて原点中心の円）
```

この**回転対称性**を使うと、$xy$ 平面全体への積分を「半径 $r$ の極細いリング（円環）への分解」に言い換えられる。

- 半径 $r$、幅 $dr$ のリングの**面積** $\approx$ （円周）× （幅）$= 2\pi r \cdot dr$
- そのリング上での**関数値**はどこも同じ：$e^{-r^2/2}$

よって：

$$\iint_{\mathbb{R}^2} e^{-(x^2+y^2)/2}\,dx\,dy = \int_0^\infty e^{-r^2/2} \cdot 2\pi r\,dr$$

これが極座標変換のヤコビアン（$dx\,dy = r\,dr\,d\theta$）の幾何学的な正体だ。$\theta$ 方向の積分（$\int_0^{2\pi}d\theta = 2\pi$）は回転対称性から自明に分離できる。

極座標変換 $x = r\cos\theta$、$y = r\sin\theta$（$x^2+y^2=r^2$、$dx\,dy = r\,dr\,d\theta$）を適用：

$$I^2 = \int_0^{2\pi}\int_0^{\infty} e^{-r^2/2}\,r\,dr\,d\theta = 2\pi \int_0^{\infty} r\,e^{-r^2/2}\,dr$$

$u = r^2/2$ と置換すると $du = r\,dr$ なので、

$$= 2\pi \int_0^{\infty} e^{-u}\,du = 2\pi\,[-e^{-u}]_0^{\infty} = 2\pi$$

よって $I = \sqrt{2\pi}$、$\int \phi(x)\,dx = I/\sqrt{2\pi} = 1$。$\square$

### 5.5 指数分布（Exponential Distribution）—— $\text{Exp}(\lambda)$

「ある出来事が起きるまでの待ち時間」の分布。$\lambda$ は単位時間あたりの平均発生回数（強度）。

$$f(x) = \lambda e^{-\lambda x} \quad (x \geq 0), \quad F(x) = 1 - e^{-\lambda x}$$

$$E[X] = \frac{1}{\lambda}, \quad \text{Var}(X) = \frac{1}{\lambda^2}$$

**特徴**：右に裾を引く歪んだ分布（正の値しかとらない）。**無記憶性**を持つ（「今まで待ったことが今後の待ち時間に影響しない」）。

**$E[X]$ の計算（部分積分）**：

$$E[X] = \int_0^{\infty} x \cdot \lambda e^{-\lambda x}\,dx = \lambda \int_0^{\infty} x e^{-\lambda x}\,dx$$

$u = x$、$dv = e^{-\lambda x}dx$（$\Rightarrow v = -e^{-\lambda x}/\lambda$）として部分積分：

$$= \lambda \left[\left[-\frac{x}{\lambda}e^{-\lambda x}\right]_0^{\infty} + \int_0^{\infty} \frac{1}{\lambda}e^{-\lambda x}\,dx\right] = \lambda \cdot \frac{1}{\lambda} \int_0^{\infty} e^{-\lambda x}\,dx = \left[-e^{-\lambda x}\right]_0^{\infty} \cdot \frac{1}{\lambda} \cdot \lambda = \frac{1}{\lambda} \quad \square$$

（$xe^{-\lambda x} \to 0$ as $x \to \infty$ を使用）

### 5.6 カイ二乗分布（Chi-squared Distribution）—— $\chi^2(k)$

$Z_1, \ldots, Z_k \overset{i.i.d.}{\sim} N(0,1)$ のとき $\chi^2 = Z_1^2 + \cdots + Z_k^2 \sim \chi^2(k)$（自由度 $k$ のカイ二乗分布）。

$$E[X] = k, \quad \text{Var}(X) = 2k$$

**計量経済学での役割**：仮説検定の検定統計量（Wald 検定、LM 検定）が漸近的にカイ二乗分布に従う（第Ⅳ部）。

### 5.7 その他の重要な分布（後の章で詳述）

| 分布 | 記号 | 主な用途 |
|-----|------|---------|
| ポアソン分布（Poisson distribution） | $\text{Poi}(\lambda)$ | 単位時間内の事象発生件数 |
| t 分布（Student's t-distribution） | $t(k)$ | 小標本での係数検定 |
| F 分布（F-distribution） | $F(k_1, k_2)$ | 分散比・重回帰の F 検定 |
| 対数正規分布（Log-normal distribution） | $\text{LN}(\mu, \sigma^2)$ | 株価・所得などの正の変数 |

---

## 6. 期待値（Expected Value）

### 6.1 定義

**定義 6.1（期待値）**：

$$E[X] = \begin{cases} \displaystyle\sum_x x \cdot p(x) & \text{（離散）} \\[8pt] \displaystyle\int_{-\infty}^{\infty} x \cdot f(x)\,dx & \text{（連続）} \end{cases}$$

ルベーグ積分で統一すると $E[X] = \int_\Omega X(\omega)\,dP(\omega)$。

> **直感**：「各値 × その値が出る確率（密度）」の合計。無限回試行したときの「長期的な平均値」の理論的な値。

**例（ベルヌーイ）**：$E[X] = 1 \cdot p + 0 \cdot (1-p) = p$

**例（標準正規）**：$xe^{-x^2/2}$ は奇関数（$f(-x) = -f(x)$）なので $\int_{-\infty}^{\infty} x \phi(x)\,dx = 0 \Rightarrow E[Z] = 0$

**例（一様分布 $U(a,b)$）**：

$$E[X] = \int_a^b x \cdot \frac{1}{b-a}\,dx = \frac{1}{b-a}\left[\frac{x^2}{2}\right]_a^b = \frac{b^2-a^2}{2(b-a)} = \frac{a+b}{2}$$

### 6.2 期待値の線形性（最重要性質）

**定理 6.1（線形性）**：

$$E[aX + bY] = a\,E[X] + b\,E[Y]$$

これは $X, Y$ が独立かどうかを問わず成立する。

*証明*（連続の場合）：

$$E[aX+bY] = \int\!\!\int (ax+by) f(x,y)\,dx\,dy = a\int\!\!\int x f(x,y)\,dx\,dy + b\int\!\!\int y f(x,y)\,dx\,dy = aE[X] + bE[Y] \quad \square$$

**関数の期待値（LOTUS）**：$Y = g(X)$ の期待値は $Y$ の分布を求めなくても、

$$E[g(X)] = \int_{-\infty}^{\infty} g(x)\,f(x)\,dx$$

で計算できる（**LOTUS = Law of the Unconscious Statistician**）。

---

## 7. 分散（Variance）

### 7.1 定義と計算公式

**定義 7.1（分散）**：

$$\text{Var}(X) = E\!\left[(X - E[X])^2\right]$$

**計算に便利な公式**：

$$\text{Var}(X) = E[X^2] - (E[X])^2$$

*展開による証明*（$\mu = E[X]$ とおく）：

$$E[(X-\mu)^2] = E[X^2 - 2\mu X + \mu^2] \overset{\text{線形性}}{=} E[X^2] - 2\mu\,E[X] + \mu^2 = E[X^2] - 2\mu^2 + \mu^2 = E[X^2] - \mu^2 \quad \square$$

> **直感**：分散は「平均からのズレの2乗」の平均。値がどれだけ平均の周りに散らばっているかを測る。2乗するのは「プラスとマイナスのズレが打ち消し合わないよう」にするため。

### 7.2 分散の性質

$$\text{Var}(aX + b) = a^2\,\text{Var}(X) \quad \text{（定数 $b$ を足しても分散は変わらない）}$$

*証明*：$E[aX+b] = aE[X]+b$ だから

$$\text{Var}(aX+b) = E[(aX+b - (aE[X]+b))^2] = E[a^2(X-E[X])^2] = a^2\,\text{Var}(X) \quad \square$$

$X \perp Y$（独立）のとき：$\text{Var}(X + Y) = \text{Var}(X) + \text{Var}(Y)$

### 7.3 各分布の期待値・分散・特徴一覧

| 分布の名前 | 記号 | $E[X]$ | $\text{Var}(X)$ | 特徴 |
|-----------|------|--------|----------------|------|
| **ベルヌーイ分布**（Bernoulli distribution） | $\text{Ber}(p)$ | $p$ | $p(1-p)$ | 0/1の2値 |
| **二項分布**（Binomial distribution） | $\text{Bin}(n,p)$ | $np$ | $np(1-p)$ | 成功回数 |
| **一様分布**（Uniform distribution） | $U(a,b)$ | $(a+b)/2$ | $(b-a)^2/12$ | 均等な確率 |
| **正規分布**（Normal distribution） | $N(\mu,\sigma^2)$ | $\mu$ | $\sigma^2$ | 釣り鐘型 |
| **指数分布**（Exponential distribution） | $\text{Exp}(\lambda)$ | $1/\lambda$ | $1/\lambda^2$ | 待ち時間 |
| **カイ二乗分布**（Chi-squared distribution） | $\chi^2(k)$ | $k$ | $2k$ | 正規の2乗和 |
| **ポアソン分布**（Poisson distribution） | $\text{Poi}(\lambda)$ | $\lambda$ | $\lambda$ | 件数（稀な事象） |

---

## 8. 共分散と相関係数

**定義 8.1（共分散）**：

$$\text{Cov}(X, Y) = E[(X - E[X])(Y - E[Y])] = E[XY] - E[X]\,E[Y]$$

**定義 8.2（相関係数）**：

$$\rho_{XY} = \frac{\text{Cov}(X, Y)}{\sqrt{\text{Var}(X)}\cdot\sqrt{\text{Var}(Y)}} \in [-1, 1]$$

$\rho > 0$：正の相関（$X$ が大きいと $Y$ も大きい傾向）  
$\rho < 0$：負の相関  
$\rho = 0$：無相関

> **重要な注意：独立と無相関は「必要十分」の関係ではない**
>
> 独立（independence）と無相関（zero correlation）の間には、**一方向のみの論理関係**しか成り立たない：
>
> $$X \perp Y \;(\text{独立}) \quad \Rightarrow \quad \text{Cov}(X,Y) = 0 \;(\text{無相関})$$
>
> しかし、**逆（$\Leftarrow$）は成り立たない**。無相関であっても独立とは限らない。
>
> | 関係 | 成立するか？ |
> |------|-------------|
> | 独立 $\Rightarrow$ 無相関 | ✅ **常に成立**（$E[XY]=E[X]E[Y]$ より） |
> | 無相関 $\Rightarrow$ 独立 | ❌ **一般に不成立** |
> | 独立 $\Leftrightarrow$ 無相関 | ❌ **必要十分ではない** |
>
> **反例**：$X \sim U(-1, 1)$（$-1$ から $1$ の一様分布）、$Y = X^2$ とする。
>
> - $Y$ は $X$ の値が決まれば完全に決まる（$X$ が分かれば $Y$ が分かる）→ **$X$ と $Y$ は完全に従属（dependent）**
> - 一方、$E[X] = 0$（対称性より）、$E[X^3] = 0$（奇関数の対称積分）だから：
>   $$\text{Cov}(X, Y) = E[XY] - E[X]\,E[Y] = E[X \cdot X^2] - 0 \cdot E[X^2] = E[X^3] = 0$$
> - したがって $\text{Cov}(X,Y) = 0$（**無相関**）だが、$X \not\perp Y$（**独立ではない**）
>
> **直感的理解**：共分散・相関係数は「**線形**な関係の強さ」しか測れない。$Y = X^2$ のような非線形な関係は捉えられないのだ。$X$ が正でも負でも $Y$ は正になるため、$X$ と $Y$ の動き方に「プラスとマイナスが対称にキャンセル」されて共分散がゼロになってしまう。

---

## 9. 確率ベクトルと分散共分散行列

### 9.1 確率ベクトル

**定義 9.1（確率ベクトル）**：$k$ 個の確率変数を縦に並べたベクトル

$$\mathbf{X} = \begin{pmatrix} X_1 \\ X_2 \\ \vdots \\ X_k \end{pmatrix}$$

を**確率ベクトル（random vector）**という。

**期待値ベクトル（平均ベクトル）**：

$$E[\mathbf{X}] = \begin{pmatrix} E[X_1] \\ \vdots \\ E[X_k] \end{pmatrix} = \boldsymbol{\mu}$$

**線形変換の期待値**：$\mathbf{Y} = A\mathbf{X} + \mathbf{b}$（$A$ は定数行列、$\mathbf{b}$ は定数ベクトル）のとき、

$$E[\mathbf{Y}] = A\,E[\mathbf{X}] + \mathbf{b} = A\boldsymbol{\mu} + \mathbf{b}$$

### 9.2 分散共分散行列

**定義 9.2（分散共分散行列）**：

$$\Sigma = \text{Var}(\mathbf{X}) = E[(\mathbf{X} - \boldsymbol{\mu})(\mathbf{X} - \boldsymbol{\mu})']$$

具体的に展開すると（$(i,j)$ 要素が $\text{Cov}(X_i, X_j)$）：

$$\Sigma = \begin{pmatrix} \text{Var}(X_1) & \text{Cov}(X_1,X_2) & \cdots & \text{Cov}(X_1,X_k) \\ \text{Cov}(X_2,X_1) & \text{Var}(X_2) & \cdots & \text{Cov}(X_2,X_k) \\ \vdots & & \ddots & \vdots \\ \text{Cov}(X_k,X_1) & \cdots & \cdots & \text{Var}(X_k) \end{pmatrix}$$

対角成分 $= $ 各変数の分散、非対角成分 $= $ 異なる変数の共分散。

**$\Sigma$ の性質**：
- 対称行列：$\Sigma = \Sigma'$（$\text{Cov}(X_i,X_j) = \text{Cov}(X_j,X_i)$ なので）
- 半正定値：$\mathbf{a}'\Sigma\mathbf{a} \geq 0$ for all $\mathbf{a}$（分散は非負のため）

**線形変換の分散公式（最重要）**：$\mathbf{Y} = A\mathbf{X} + \mathbf{b}$ のとき

$$\text{Var}(\mathbf{Y}) = A\,\Sigma\,A'$$

*証明*：$\boldsymbol{\mu}_Y = A\boldsymbol{\mu} + \mathbf{b}$ なので、

$$\text{Var}(\mathbf{Y}) = E[(\mathbf{Y}-\boldsymbol{\mu}_Y)(\mathbf{Y}-\boldsymbol{\mu}_Y)'] = E[(A(\mathbf{X}-\boldsymbol{\mu}))(A(\mathbf{X}-\boldsymbol{\mu}))'] = A\,E[(\mathbf{X}-\boldsymbol{\mu})(\mathbf{X}-\boldsymbol{\mu})']\,A' = A\Sigma A' \quad \square$$

> **第7章への接続**：OLS 推定量は $\hat{\beta} = (X'X)^{-1}X'y$ と書ける。これは $y = X\beta + \varepsilon$ に代入すると $\hat{\beta} = \beta + (X'X)^{-1}X'\varepsilon$ となる。$A = (X'X)^{-1}X'$、$\mathbf{X} = \varepsilon$ として上の公式を使うと $\text{Var}(\hat{\beta}) = \sigma^2(X'X)^{-1}$ が導かれる。

---

## 10. 同時分布・周辺分布・条件付き分布

### 10.1 同時累積分布関数（Joint CDF）

**定義 10.1（同時 CDF）**：

$$F(x, y) = P(X \leq x \text{ かつ } Y \leq y)$$

**性質**：$F(-\infty, y) = F(x, -\infty) = 0$、$F(\infty, \infty) = 1$

### 10.2 同時確率密度関数（Joint PDF）

**定義 10.2（同時 PDF）**：関数 $f(x, y) \geq 0$ が

$$P(a \leq X \leq b,\ c \leq Y \leq d) = \int_a^b\!\!\int_c^d f(x, y)\,dy\,dx$$

を満たすとき、$f$ を $(X, Y)$ の**同時確率密度関数（Joint PDF）**という。

規格化条件：$\int_{-\infty}^{\infty}\!\int_{-\infty}^{\infty} f(x,y)\,dy\,dx = 1$

同時 CDF との関係：$f(x,y) = \dfrac{\partial^2 F(x,y)}{\partial x\,\partial y}$

### 10.3 周辺分布（Marginal Distribution）

同時分布から片方の変数を「積分で消す」と周辺分布が得られる。

**定義 10.3（周辺 PDF）**：

$$f_X(x) = \int_{-\infty}^{\infty} f(x,y)\,dy \quad \text{（$y$ を積分して消す）}$$

$$f_Y(y) = \int_{-\infty}^{\infty} f(x,y)\,dx \quad \text{（$x$ を積分して消す）}$$

```
【周辺分布のイメージ】

  y
  │ 同時分布 f(x,y) の
  │ 「山」を想像する。
  │
  │     ／＼
  │    /   \
  └──────────── x
     ↓↓↓（y方向に押し潰す）
  f_X(x)：x方向の輪郭（周辺分布）
```

### 10.4 条件付き分布（Conditional Distribution）

**定義 10.4（条件付き PDF）**：$f_X(x) > 0$ のとき、

$$f_{Y \mid X}(y \mid x) = \frac{f(x, y)}{f_X(x)}$$

これは $y$ を変数としたとき正しい PDF になっている：$f_{Y|X}(y|x) \geq 0$、$\int f_{Y|X}(y|x)\,dy = 1$

**乗法定理（密度版）**：

$$f(x, y) = f_{Y \mid X}(y \mid x) \cdot f_X(x)$$

**独立性の密度版**：

$$X \perp Y \iff f(x, y) = f_X(x) \cdot f_Y(y) \quad \text{（全ての $x, y$ で成立）}$$

---

## 11. 条件付き期待値（Conditional Expectation）

### 11.1 定義

**定義 11.1（条件付き期待値）**：

$$E[Y \mid X = x] = \int_{-\infty}^{\infty} y \cdot f_{Y \mid X}(y \mid x)\,dy \quad \text{（連続）}$$

**重要な区別**：

- $E[Y \mid X = x]$：$X$ が具体的な値 $x$ をとるという条件のもとでの $Y$ の期待値。これは $x$ の**関数（数値を返す）**。
- $E[Y \mid X]$：確率変数 $X$ を与えたときの条件付き期待値。これは $X$ の**確率変数**。

**例**：$Y = 2X + \varepsilon$（$\varepsilon \sim N(0,1)$、$\varepsilon \perp X$）のとき、

$$E[Y \mid X = x] = 2x + E[\varepsilon \mid X = x] = 2x + 0 = 2x$$

### 11.2 重要な性質

**CE1（繰り返し期待値の法則・LIE）**：

$$E\!\left[E[Y \mid X]\right] = E[Y]$$

LIE = Law of Iterated Expectations（繰り返し期待値の法則）。

*証明*（連続）：

$$E[E[Y|X]] = \int_{-\infty}^{\infty} E[Y|X=x]\cdot f_X(x)\,dx = \int\!\!\int y\cdot f_{Y|X}(y|x)\,f_X(x)\,dy\,dx = \int\!\!\int y\cdot f(x,y)\,dy\,dx = E[Y] \quad \square$$

> **例**：男性の平均賃金が 50万円、女性が 40万円で、男女比が 6:4 とする。全体の平均賃金 = $50 \times 0.6 + 40 \times 0.4 = 46$ 万円。これが LIE の意味。

**CE2（予測不可能性）**：$E[\varepsilon \mid X] = 0 \Rightarrow E[X\varepsilon] = 0$

*証明*：LIE を使って

$$E[X\varepsilon] = E\!\left[E[X\varepsilon \mid X]\right] = E\!\left[X \cdot E[\varepsilon \mid X]\right] = E[X \cdot 0] = 0 \quad \square$$

これは OLS の仮定「$E[\varepsilon \mid X] = 0$（CLM3）」の数学的帰結であり、$\varepsilon$ と $X$ が単なる「無相関」を超えた「直交性」を持つことを意味する。

**CE3（線形性）**：$E[aY + bZ \mid X] = a\,E[Y \mid X] + b\,E[Z \mid X]$

**CE4（既知変数の引き出し）**：$h(X)$ が $X$ の関数ならば

$$E[h(X) \cdot Y \mid X] = h(X) \cdot E[Y \mid X]$$

---

## 12. 独立同分布（i.i.d.）

**定義 12.1（独立同分布・i.i.d.）**：

確率変数の列 $X_1, X_2, \ldots, X_n$ が **i.i.d.（independent and identically distributed）**、すなわち「独立同分布（どくりつどうぶんぷ）」であるとは、

1. **相互独立**：任意の部分集合 $\{X_{i_1}, \ldots, X_{i_k}\}$ の結合分布が周辺分布の積に等しい
2. **同分布**：全ての $i$ について $X_i$ が同じ分布 $F$ に従う

の両方が成立することをいう。表記：$X_1, \ldots, X_n \overset{i.i.d.}{\sim} F$。

> **直感**：「全く同じコインを $n$ 回投げる」。各回の結果は前の結果に影響されず（独立）、どの回も同じコインを使う（同分布）。

**i.i.d. が崩れる状況と対処**：

| 状況 | 独立性が崩れる理由 | 対処法（本書での登場） |
|-----|-----------------|-------------------|
| 時系列データ | $X_t$ と $X_{t-1}$ が相関（自己相関） | HAC 標準誤差（第Ⅲ部） |
| パネルデータ | 同一個体の観測が相関（クラスター相関） | クラスター標準誤差（第Ⅲ部） |
| 不均一分散 | $\text{Var}(X_i)$ が $i$ によって異なる | HC 標準誤差（第Ⅲ部） |

---

## 13. 大数の法則（Law of Large Numbers）

### 13.1 問題設定

$X_1, X_2, \ldots$ が i.i.d. で $E[X_i] = \mu$ とする。標本平均

$$\bar{X}_n = \frac{1}{n}\sum_{i=1}^n X_i$$

は $n$ が増えるとどうなるか？直感的には「真の平均 $\mu$ に近づく」はずだが、それを数学的に示すのが大数の法則だ。

### 13.2 マルコフの不等式（Markov's Inequality）

大数の法則の証明の出発点となる不等式。

**定理 13.1（マルコフの不等式）**：$Y \geq 0$ の確率変数で $E[Y] < \infty$ のとき、任意の $c > 0$ に対して、

$$P(Y \geq c) \leq \frac{E[Y]}{c}$$

**直感的な理解**：

```
  f_Y(y)
   │                          ← Y の密度関数
   │    ███████████████
   │    █████████████████
   │    ████████████████████
   └──────────┼─────────── y
              c

  Y ≥ c の部分（右側）の確率は、
  E[Y] = ∫y f(y)dy の中で
  「y ≥ c の部分」だけを取っても E[Y] 以下のはず。
```

*証明*：$Y \geq 0$ だから

$$E[Y] = \int_0^{\infty} y\,f_Y(y)\,dy \geq \int_c^{\infty} y\,f_Y(y)\,dy \geq \int_c^{\infty} c\,f_Y(y)\,dy = c \cdot P(Y \geq c)$$

両辺を $c > 0$ で割って $P(Y \geq c) \leq E[Y]/c$。$\square$

**具体例**：ある学校の生徒の平均点が 60点（$E[Y] = 60$）とする。「90点以上の生徒は何割か？」

$$P(Y \geq 90) \leq \frac{60}{90} = \frac{2}{3} \approx 66.7\%$$

これは「上界」（本当の確率はこれ以下）を与えるだけで、緩い不等式である。しかし仮定が少ない（$Y \geq 0$ と $E[Y] < \infty$ だけ）のが強み。

### 13.3 チェビシェフの不等式（Chebyshev's Inequality）

マルコフの不等式を「分散を使って改良」した不等式。

**定理 13.2（チェビシェフの不等式）**：$E[X] = \mu$, $\text{Var}(X) = \sigma^2$ のとき、任意の $\varepsilon > 0$ に対して、

$$P(|X - \mu| \geq \varepsilon) \leq \frac{\sigma^2}{\varepsilon^2}$$

同値的に：$P(|X - \mu| < \varepsilon) \geq 1 - \frac{\sigma^2}{\varepsilon^2}$

**直感**：分散（散らばりの大きさ）が小さいほど、平均から大きく外れる確率は小さい。

*証明*：$Y = (X - \mu)^2$ とおく。$Y \geq 0$ かつ $E[Y] = \sigma^2$。$c = \varepsilon^2$ としてマルコフの不等式を適用：

$$P\!\left((X-\mu)^2 \geq \varepsilon^2\right) \leq \frac{E[(X-\mu)^2]}{\varepsilon^2} = \frac{\sigma^2}{\varepsilon^2}$$

$(X-\mu)^2 \geq \varepsilon^2 \iff |X-\mu| \geq \varepsilon$ だから結論が得られる。$\square$

**具体例で感覚をつかむ**：

$\mu = 0$、$\sigma^2 = 1$（例：標準正規分布に限らず、平均0・分散1の任意の分布）。

| $\varepsilon$（平均からの距離） | 上界 $\sigma^2/\varepsilon^2$ | 正規分布の実際の確率 |
|---|---|---|
| $\varepsilon = 1$ | $1/1 = 100\%$ | $\approx 31.7\%$ |
| $\varepsilon = 2$ | $1/4 = 25\%$ | $\approx 4.6\%$ |
| $\varepsilon = 3$ | $1/9 \approx 11.1\%$ | $\approx 0.3\%$ |

チェビシェフの上界は緩い（実際の確率よりずっと大きい）が、「分布の形を問わずに使える」という汎用性がある。

### 13.4 弱大数の法則（WLLN: Weak Law of Large Numbers）

**定理 13.3（弱大数の法則）**：$X_i \overset{i.i.d.}{\sim}$、$E[X_i] = \mu$、$\text{Var}(X_i) = \sigma^2 < \infty$ のとき、任意の $\varepsilon > 0$ に対して

$$P\!\left(|\bar{X}_n - \mu| \geq \varepsilon\right) \to 0 \quad (n \to \infty)$$

これを「$\bar{X}_n$ が $\mu$ に**確率収束（convergence in probability）**する」といい、$\bar{X}_n \overset{p}{\to} \mu$ と書く。

*証明*：

まず $\bar{X}_n$ の期待値と分散を計算する。

$$E[\bar{X}_n] = E\!\left[\frac{1}{n}\sum_{i=1}^n X_i\right] = \frac{1}{n}\sum_{i=1}^n E[X_i] = \frac{1}{n} \cdot n\mu = \mu$$

$$\text{Var}(\bar{X}_n) = \text{Var}\!\left(\frac{1}{n}\sum_{i=1}^n X_i\right) \overset{\text{独立性}}{=} \frac{1}{n^2}\sum_{i=1}^n \text{Var}(X_i) = \frac{1}{n^2} \cdot n\sigma^2 = \frac{\sigma^2}{n}$$

チェビシェフの不等式を $\bar{X}_n$（平均 $\mu$、分散 $\sigma^2/n$）に適用：

$$P(|\bar{X}_n - \mu| \geq \varepsilon) \leq \frac{\sigma^2/n}{\varepsilon^2} = \frac{\sigma^2}{n\varepsilon^2}$$

$n \to \infty$ のとき右辺 $\to 0$（$\sigma^2, \varepsilon$ は固定）。よって $P(|\bar{X}_n - \mu| \geq \varepsilon) \to 0$。$\square$

### 13.5 強大数の法則（SLLN: Strong Law of Large Numbers）

**定理 13.4（強大数の法則）**：$X_i \overset{i.i.d.}{\sim}$、$E[|X_i|] < \infty$ のとき（$E[X_i] = \mu$）、

$$P\!\left(\lim_{n \to \infty} \bar{X}_n = \mu\right) = 1$$

これを「$\bar{X}_n$ が $\mu$ に**概収束（almost sure convergence）**する」といい、$\bar{X}_n \overset{a.s.}{\to} \mu$ と書く。

**WLLN と SLLN の違い**：

- WLLN（弱）：「各 $n$ を固定したとき、$\bar{X}_n$ が $\mu$ から離れる確率が 0 に近づく」
- SLLN（強）：「$n$ を無限に増やしていくとき、$\bar{X}_n$ の列そのものが $\mu$ に収束する確率が 1」

SLLN の方が強い主張で、証明も難しい（WLLN は分散の有限性が必要だが、SLLN は期待値の有限性だけで十分）。

---

## 14. 中心極限定理（Central Limit Theorem）

### 14.1 動機

大数の法則は「$\bar{X}_n$ がどの値（$\mu$）に収束するか」を教えた。では「$\bar{X}_n$ は $\mu$ の周りにどんな形で散らばっているか」——つまり $\bar{X}_n$ の**分布の形**はどうなっているか？

CLT（中心極限定理）の答え：「$n$ が十分大きければ、$\bar{X}_n$ の（標準化した）分布は正規分布に近づく」。元の $X_i$ の分布が何であれ。

### 14.2 中心極限定理（CLT）

**定理 14.1（中心極限定理・CLT）**：$X_i \overset{i.i.d.}{\sim}$、$E[X_i] = \mu$、$\text{Var}(X_i) = \sigma^2 < \infty$ のとき、

$$Z_n = \frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} = \frac{\sqrt{n}(\bar{X}_n - \mu)}{\sigma} \overset{d}{\to} N(0, 1) \quad (n \to \infty)$$

$\overset{d}{\to}$ は**分布収束（convergence in distribution）**：任意の $z$ に対して $P(Z_n \leq z) \to \Phi(z)$。

**使いやすい形**：$n$ が十分大きいとき近似的に

$$\bar{X}_n \overset{\text{近似}}{\sim} N\!\left(\mu,\ \frac{\sigma^2}{n}\right), \qquad \sum_{i=1}^n X_i \overset{\text{近似}}{\sim} N(n\mu,\ n\sigma^2)$$

> **驚くべき点**：$X_i$ の分布が一様分布でも指数分布でもベルヌーイ分布でも、$n$ が大きければ標本平均の分布は正規分布に近づく。これが「正規分布があらゆる場所に現れる」理由の正体だ。

### 14.3 CLT の証明（モーメント母関数法）

モーメント母関数（Moment Generating Function, MGF）を使った証明を示す。

**MGF の定義**：確率変数 $X$ の MGF を

$$M_X(t) = E[e^{tX}]$$

と定義する。MGF が重要な理由は「分布が MGF によって一意に決まる」からだ（同じ MGF を持つ分布は同じ）。

**ステップ 1：標準化**

$Y_i = \dfrac{X_i - \mu}{\sigma}$ とおく（$E[Y_i] = 0$, $\text{Var}(Y_i) = 1$）。

$$Z_n = \frac{1}{\sqrt{n}}\sum_{i=1}^n Y_i$$

**ステップ 2：$Y_i$ の MGF をマクローリン展開する**

$M_Y(t) = E[e^{tY}]$ をマクローリン展開（0.3節）する。

$e^{tY} = 1 + tY + \dfrac{(tY)^2}{2!} + \dfrac{(tY)^3}{3!} + \cdots$ の両辺の期待値をとると：

$$M_Y(t) = E[e^{tY}] = 1 + t\,E[Y] + \frac{t^2}{2}E[Y^2] + \frac{t^3}{6}E[Y^3] + \cdots$$

$E[Y] = 0$、$E[Y^2] = \text{Var}(Y) = 1$ を代入：

$$M_Y(t) = 1 + \frac{t^2}{2} + \frac{t^3}{6}E[Y^3] + \cdots = 1 + \frac{t^2}{2} + O(t^3)$$

（$O(t^3)$ は「$t^3$ オーダーの残余項」を表す記号）

**ステップ 3：$Z_n$ の MGF を計算する**

独立性（$Y_i$ 同士が独立）より、独立な確率変数の和の MGF は MGF の積：

$$M_{Z_n}(t) = E\!\left[e^{tZ_n}\right] = E\!\left[e^{t \cdot \frac{1}{\sqrt{n}}\sum Y_i}\right] = \prod_{i=1}^n E\!\left[e^{\frac{t}{\sqrt{n}}Y_i}\right] = \left(M_Y\!\left(\frac{t}{\sqrt{n}}\right)\right)^n$$

**ステップ 4：$M_Y(t/\sqrt{n})$ を展開する**

ステップ 2 の結果で $t$ を $t/\sqrt{n}$ に置き換える：

$$M_Y\!\left(\frac{t}{\sqrt{n}}\right) = 1 + \frac{1}{2}\left(\frac{t}{\sqrt{n}}\right)^2 + O\!\left(n^{-3/2}\right) = 1 + \frac{t^2}{2n} + O\!\left(n^{-3/2}\right)$$

**ステップ 5：$n \to \infty$ の極限をとる**

$$M_{Z_n}(t) = \left(1 + \frac{t^2}{2n} + O\!\left(n^{-3/2}\right)\right)^n$$

0.3節の極限公式 $\left(1 + \dfrac{c}{n}\right)^n \to e^c$ で $c = t^2/2$ とすると：

$$M_{Z_n}(t) \to e^{t^2/2} \quad (n \to \infty)$$

**ステップ 6：$N(0,1)$ の MGF と一致することを確認**

$Z \sim N(0,1)$ の MGF：

$$M_Z(t) = E[e^{tZ}] = \int_{-\infty}^{\infty} e^{tx} \cdot \frac{1}{\sqrt{2\pi}}e^{-x^2/2}\,dx = \frac{1}{\sqrt{2\pi}}\int e^{-(x^2-2tx)/2}\,dx$$

指数部分を平方完成：$x^2 - 2tx = (x-t)^2 - t^2$

$$= \frac{e^{t^2/2}}{\sqrt{2\pi}}\int_{-\infty}^{\infty} e^{-(x-t)^2/2}\,dx = e^{t^2/2} \cdot 1 = e^{t^2/2}$$

$Z_n$ の MGF $\to e^{t^2/2}$ = $N(0,1)$ の MGF。MGF が収束すれば分布も収束するので、$Z_n \overset{d}{\to} N(0,1)$。$\square$

### 14.4 ベリー＝エッセンの定理（Berry-Esseen Theorem）

CLT は「$n \to \infty$ での近似」を保証するが、有限の $n$ での近似の誤差はどのくらいか？これに答えるのが**ベリー＝エッセンの定理**だ。

**定理 14.2（ベリー＝エッセンの定理）**：$X_i \overset{i.i.d.}{\sim}$、$E[X_i] = \mu$、$\text{Var}(X_i) = \sigma^2$、かつ $E[|X_i - \mu|^3] = \rho < \infty$ のとき、全ての実数 $z$ と全ての $n$ に対して、

$$\sup_{z \in \mathbb{R}} \left| P(Z_n \leq z) - \Phi(z) \right| \leq \frac{C \cdot \rho}{\sigma^3 \sqrt{n}}$$

ここで $C$ は普遍定数（どんな分布でも共通）で、$C \leq 0.4748$（最良の上界として知られている）。

**定理の意味を読み解く**：

- $\sup_{z}|P(Z_n \leq z) - \Phi(z)|$：CLT の近似誤差の「最大値」。どんな $z$ の点でも近似誤差はこの値以下。
- 誤差は $O(1/\sqrt{n})$ のオーダーで減少する——$n$ を4倍にすると誤差は半分になる。
- $\rho/\sigma^3$ は「3次モーメントと分散の比」。元の分布が正規分布に近いほど（歪みが小さいほど）小さく、収束が速い。

**具体的な計算例**：

ベルヌーイ分布 $\text{Ber}(0.5)$（コイン投げ）の場合：
- $\mu = 0.5$, $\sigma^2 = 0.25$, $E[|X-\mu|^3] = 0.5 \times |0.5|^3 + 0.5 \times |{-0.5}|^3 = 0.125 = \rho$

$$\text{誤差} \leq \frac{0.4748 \times 0.125}{0.5^3 \times \sqrt{n}} = \frac{0.05935}{0.125\sqrt{n}} \approx \frac{0.4748}{\sqrt{n}}$$

| $n$ | 誤差の上界 | CLT の近似精度 |
|-----|-----------|-------------|
| $n = 10$ | $\approx 0.150$ | 最大15%の誤差あり |
| $n = 30$ | $\approx 0.087$ | 最大8.7%の誤差あり |
| $n = 100$ | $\approx 0.047$ | 最大4.7%の誤差あり |
| $n = 1000$ | $\approx 0.015$ | 最大1.5%の誤差あり |

**実践的な教訓**：

1. 「$n \geq 30$ なら CLT が使える」という経験則は目安に過ぎない。元の分布が正規に近いなら $n = 5$ でも十分なことがある。元の分布が極端に歪んでいれば $n = 1000$ でも不十分なことがある。

2. 誤差は$O(1/\sqrt{n})$ なので、精度を10倍上げたければ $n$ を100倍必要がある。

3. 計量経済学での実践：推定量が正規分布に近いかどうかは、理論的な議論に加えてシミュレーション（モンテカルロ法）で確認するのがよい。

---

## 15. 批判層：仮定の限界と計量経済学への接続

### 15.1 確率の2大解釈

| | 頻度論（Frequentist） | ベイズ（Bayesian） |
|---|---|---|
| 確率の意味 | 長期的な相対頻度 | 信念の度合い（degree of belief） |
| パラメータの扱い | 真の固定値（確率変数ではない） | 確率変数（事前分布を持つ） |
| 推論の手段 | OLS・MLE・信頼区間・p値 | 事後分布・信用区間（credible interval） |

本書は第Ⅸ部まで主に頻度論に基づく。

### 15.2 本章で学んだ概念と計量経済学の対応

| 本章の概念 | 計量経済学での役割 | 登場する章 |
|-----------|----------------|----------|
| $E[\varepsilon \mid X] = 0$ | OLS 不偏性の前提（CLM3） | 第7章 |
| $\text{Var}(\varepsilon \mid X) = \sigma^2$ | 等分散仮定（CLM5） | 第7・8章 |
| $\Sigma = E[\mathbf{X}\mathbf{X}'] $（分散共分散行列） | $(X'X)^{-1}$ の漸近的収束先 | 第11章 |
| LIE（繰り返し期待値の法則） | $E[\hat{\beta}] = \beta$ の証明 | 第7章 |
| WLLN | $\hat{\beta} \overset{p}{\to} \beta$（一致性） | 第11章 |
| CLT | $\sqrt{n}(\hat{\beta}-\beta) \overset{d}{\to} N(0, \sigma^2 Q^{-1})$ | 第11章 |
| ベリー＝エッセン定理 | 有限標本での検定の精度評価 | 第Ⅳ部 |

---

## 16. 実装（Python / R / Stata）

### 16.1 Python

```python
import numpy as np
import matplotlib.pyplot as plt
import japanize_matplotlib
from scipy import stats

rng = np.random.default_rng(seed=42)

# -------------------------------------------------------
# 1. 大数の法則：様々な分布で収束を確認
# -------------------------------------------------------
n_max = 5000
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, (name, samples, true_mean) in zip(axes, [
    ("一様分布 U(0,1)",   rng.uniform(0, 1, n_max),           0.5),
    ("指数分布 Exp(1)",   rng.exponential(1, n_max),           1.0),
    ("ベルヌーイ p=0.3",  (rng.uniform(0,1,n_max) < 0.3)*1.0, 0.3),
]):
    cumavg = np.cumsum(samples) / np.arange(1, n_max + 1)
    ax.plot(cumavg, lw=1.2, color="steelblue")
    ax.axhline(true_mean, color="red", lw=1.5, linestyle="--", label=f"真の平均 μ={true_mean}")
    ax.set_title(name); ax.set_xlabel("試行回数 n（対数スケール）")
    ax.set_ylabel("標本平均"); ax.set_xscale("log"); ax.legend(fontsize=9)
plt.suptitle("大数の法則：どの分布でも標本平均は真の期待値に収束する")
plt.tight_layout(); plt.savefig("lln.png", dpi=150); plt.show()

# -------------------------------------------------------
# 2. 中心極限定理：指数分布（歪んだ分布）でも n が増えれば正規に近づく
# -------------------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
x_grid = np.linspace(-4, 4, 300)

for ax, n in zip(axes, [1, 5, 30, 100]):
    means = np.array([rng.exponential(1, n).mean() for _ in range(10000)])
    z = (means - 1) / (1 / np.sqrt(n))  # 標準化
    ax.hist(z, bins=60, density=True, alpha=0.6, color="steelblue", label=f"n={n}")
    ax.plot(x_grid, stats.norm.pdf(x_grid), "r-", lw=2, label="N(0,1)")
    ax.set_title(f"n = {n}"); ax.legend(fontsize=8); ax.set_xlim(-4, 4)
plt.suptitle("CLT：指数分布（右に歪んだ分布）の標本平均の分布")
plt.tight_layout(); plt.savefig("clt.png", dpi=150); plt.show()

# -------------------------------------------------------
# 3. ベリー＝エッセン：n と近似誤差の関係
# -------------------------------------------------------
ns = np.array([10, 30, 100, 300, 1000, 3000])
# ベルヌーイ p=0.5：ρ/σ³ = 0.125/0.125 = 1
upper_bounds = 0.4748 / np.sqrt(ns)

plt.figure(figsize=(7, 4))
plt.loglog(ns, upper_bounds, "b-o", label="ベリー＝エッセン上界")
plt.xlabel("標本サイズ n"); plt.ylabel("近似誤差の上界")
plt.title("ベリー＝エッセン定理：n が増えると誤差は 1/√n で減少")
plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout()
plt.savefig("berry_esseen.png", dpi=150); plt.show()

# -------------------------------------------------------
# 4. チェビシェフの不等式の可視化
# -------------------------------------------------------
eps_vals = np.linspace(0.1, 5, 200)
sigma = 1.0
chebyshev = np.minimum(sigma**2 / eps_vals**2, 1.0)
normal_prob = 2 * stats.norm.sf(eps_vals / sigma)

plt.figure(figsize=(7, 4))
plt.plot(eps_vals, chebyshev, "b-", lw=2, label="チェビシェフ上界 σ²/ε²")
plt.plot(eps_vals, normal_prob, "r--", lw=2, label="N(0,1) の実際の確率")
plt.xlabel("ε（平均からの距離）"); plt.ylabel("P(|X-μ| ≥ ε)")
plt.title("チェビシェフの不等式：上界（青）は実際の確率（赤）より大きい")
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig("chebyshev.png", dpi=150); plt.show()
```

### 16.2 R

```r
set.seed(42)

# 大数の法則（指数分布）
n_max <- 5000
x <- rexp(n_max, rate=1)
cum_avg <- cumsum(x) / seq_along(x)
plot(cum_avg, type="l", log="x", col="steelblue", lwd=1.5,
     main="大数の法則（指数分布）", xlab="n（対数）", ylab="標本平均")
abline(h=1, col="red", lwd=2, lty=2)

# CLT（指数分布）
par(mfrow=c(1,4))
for (n in c(1, 5, 30, 100)) {
  means <- replicate(10000, mean(rexp(n, 1)))
  z <- (means - 1) / (1/sqrt(n))
  hist(z, prob=TRUE, breaks=50, col="lightblue", main=paste("n =", n),
       xlab="標準化標本平均", xlim=c(-4,4))
  curve(dnorm, from=-4, to=4, col="red", lwd=2, add=TRUE)
}

# ベリー＝エッセン上界
ns <- c(10, 30, 100, 300, 1000, 3000)
upper <- 0.4748 / sqrt(ns)
plot(ns, upper, log="xy", type="b", pch=16,
     main="ベリー＝エッセン：誤差上界", xlab="n", ylab="上界")
```

### 16.3 Stata

```stata
* 大数の法則
clear; set obs 5000; set seed 42
gen x = rexponential(1)
gen cum_mean = sum(x) / _n; gen n = _n
twoway line cum_mean n, xscale(log) yline(1, lp(dash) lc(red)) ///
  title("大数の法則（指数分布）") xtitle("n（対数）") ytitle("標本平均")

* CLT（n=30 での確認）
clear; set obs 10000; set seed 42
forvalues i = 1/30 { gen x`i' = rexponential(1) }
egen xbar = rowmean(x*)
gen z = (xbar - 1) / (1/sqrt(30))
histogram z, normal bin(50) xtitle("標準化標本平均") ///
  title("CLT：指数分布、n=30")
```

---

## 17. 章末問題

### レベル1（確認）

1. 「標本空間」「実数直線」「確率測度」をそれぞれ言葉で説明し、具体例（さいころ / 身長の測定）を使って対応関係を示せ。

2. 「確率密度関数（PDF）」と「確率質量関数（PMF）」の違いを説明せよ。特に「$f(x) = 3$ でも構わない」のはなぜか。

3. $f(x) = 2x$（$0 \leq x \leq 1$）が正しい PDF かどうか確認し、$P(0.5 \leq X \leq 1)$ を求めよ。

4. $X \sim U(0, 4)$ のとき $E[X]$、$\text{Var}(X)$、$P(1 \leq X \leq 3)$、$F(2)$ を求めよ。

5. $Z \sim N(0,1)$ のとき $P(-1.96 \leq Z \leq 1.96)$ を求めよ（数表を使う）。

6. 分散の公式 $\text{Var}(X) = E[X^2] - (E[X])^2$ を展開から導け。

7. LIE（繰り返し期待値の法則）を言葉で説明し、「男女別の平均賃金と全体の平均賃金」の例で説明せよ。

8. マクローリン展開で $e^x$ を 3 次まで展開せよ。

### レベル2（応用）

1. $X_1, \ldots, X_n \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$ のとき $\bar{X}_n$ の分布を導け。（ヒント：正規分布の再生性）

2. $X \sim \text{Exp}(\lambda)$ について $E[X]$ と $E[X^2]$ を部分積分で計算し、$\text{Var}(X) = 1/\lambda^2$ を示せ。

3. $Y = 2X + 3$ のとき $E[Y]$、$\text{Var}(Y)$ を $E[X]$、$\text{Var}(X)$ で表せ。

4. 同時密度 $f(x,y) = 6xy$（$0 \leq x \leq 1$, $0 \leq y \leq \sqrt{1-x}$、他は 0）の周辺密度 $f_X(x)$ を求めよ。

5. $X \perp Y$（独立）ならば $\text{Cov}(X,Y) = 0$ を示せ。逆が成り立たない反例を構成せよ。

6. チェビシェフの不等式を使って WLLN を証明せよ（$\sigma^2 < \infty$ を仮定）。

### レベル3（批判）

1. 「$n \geq 30$ あれば CLT が使える」という主張をベリー＝エッセンの定理を用いて批判的に検討せよ。どのような分布では $n = 30$ では不十分か。

2. 以下の論文記述を批判せよ：「誤差項の相関係数は 0.03 であり統計的に非有意（p = 0.12）だった。よって誤差項は独立と判断し OLS を実施した。」（ヒント：$\text{Cov}=0$ と独立の違い、帰無仮説の棄却不能の意味）

3. 「$E[\varepsilon \mid X] = 0$ は $E[\varepsilon] = 0$ より強い仮定である」ことを証明し、なぜ計量経済学では前者が必要なのかを論じよ。（ヒント：LIE を使え）

4. マクローリン展開を用いて $(1 + x/n)^n \to e^x$（$n \to \infty$）が成り立つことを示せ。（ヒント：両辺の対数をとる）

---

## 章末注

**確率の公理的定義**（コルモゴロフ, 1933）は現代確率論・統計学の基盤。

**ルベーグ積分**の完全な構成（単純関数による近似→非負関数→一般関数）は付録Cで行う。

**大数の法則と CLT の関係**：LLN は「平均が収束する先（$\mu$）」を教え、CLT は「収束の速さと形状（$\sim N(0,1)$）」を教える。2つがそろって初めて「点推定」と「区間推定」が成立する。

**ベリー＝エッセンの定数 $C$**：最良の上界として $C \leq 0.4748$ が 2010年代に確立された（Shevtsova, 2010）。理論的な下界は $C \geq (\sqrt{10}+3)/(6\sqrt{2\pi}) \approx 0.4097$ で、真の値はこの範囲にある。

---

*← 前の章：なし（第1章）*
*→ 次の章：[第2章：確率変数と分布（発展）](ch02_確率変数と分布.md)*
