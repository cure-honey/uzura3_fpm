# uzura3_fpm

MP3 encoder in Fortran95.

command line 引数を取る所と stream access 以外は、Fortran95 で書いてあると思います。（これらは Fortran 2003 になります。)

15 年放置の後、形式的エンコーダー部分のバグが取れました。その結果、心理音響部もバグが取ることが出来ました。

思いのほか良い音が出るようになり、mp3 128 kbps で CD 音質を出せるという主張をようやく信じる気になりましたｗ

## bug 

mp3 のプログラムは scale dactor を探すところが高次ベクトル空間の残差ノルム最小化の問題になっているのですが、時間がかかるので早めに探索を諦める if 文を１個入れていたのですが、それが筋が悪くて悪さをしていました。最初 Cerelon 300A 機で作り始めて実行がとても遅かったので、思い付きで入れてしまった･･･

これを外した後は、理屈通りにパラメータ変化に従って動くようになり、音を聞くことで心理音響解析部のバグも修正できました。

心理音響部については、そもそも dB とか音響学の理解が全くできておらず Power と音圧を混同して単位系が狂っているなど、プログラム以前に内容が間違っていました。ただ基本的なプログラム構造は問題ありませんでしたのでサブルーチン中の修正で済んだのが不幸中の幸いでした。

以前のプログラムは恥ずかしいので消し去りたいところですｗ

## Psychoacoustic model  

心理音響解析のモデルは Bosse Lincoln の "An Experimental High Fidelity Perceptual Audio Coder Project in MUS420 Win 97" (1998) に従っています。ただマスキングの効果は線形和ではなく、 mp1/mp2 encoder に採用したもののように非線形和を取っています。論文中でも非線形和について触れられていますが、モデルとしては採用されていません。非線形和の指数は Lincoln に従って 0.5 にしてあります。（成分のルートの和を取って二乗する 1/2-norm に相当します。） 

なお scale factor を回すよりも 1 frame 中の granule/channel 間の bit 分配がとても重要になります。上記論文は独自のフォーマットを使っているので、この部分は uzura 独自のものです。かなり適当に作りましたが極めて重要で、もっとまじめにやればさらに音質は向上すると思います。

## mixed block 

デフォルトでは long block, short block の他に、これらの二つの折衷で、低周波バンドは long block 高周波バンドは short block を利用する mixed block も使うようになっています。この点が特徴と言えるでしょうか。しかしこれは音質にさほど寄与しません。

ところで mixed block に対応していない再生機器やソフトがあります、この時 mixed block に切り替わる所でガザガザという音が少しします。そういう場合には mixed block を使わず long block と short block しか使わないようにするオプション -s を利用してください。

SONY Walkman NW A17 は対応、SONY Xperia 8 は非対応でした。

## 実行例 fortran package manager (fpm) 利用例

実行は Fortran Package Manager (fpm) を利用すると簡単にできます。

心理音響解析の非線形和に時間がかかるので、Haswell CPU 3.7GHz だと、エンコードするのに曲 1 分当たり 3 分くらいかかります。論文の式を忠実に実装して、最適化の類はしていません。心理音響解析部は精度もいらないし計算ヘビーでメモリーアクセスは余りないので高速化の余地は大いにあると思います。

** file_name.wa --> file_name.mp3 **

- debug run

 fpm run -- file_name

- release run

 fpm run --profile release -- file_name

- run with options: bitrate 320kBps, normal stereo, no mix-block s

 fpm run --profile release -- -b 14 -ns -s file_name

- non-linear sum; p-norm for adding spreading function  

 fpm run --profile release -- -pow 0.4 -- file_name 

- lineaar sum; 1-norm

 fpm run --profile release -- -pow 1.0 -- file_name 

- lineaar sum; 1-norm (another formulation)

 fpm run --profile release -- -1norm -- file_name 

### Mid-Side stereo

 normal stereo と mid-side の切り替えがうまく行かない時があります。オプション xms 0.0 で normal stereo，xms 1.0 で MS になりやすくなります

- normal stereo 多め

 fpm run --profile release -- -xms 0.2  file_name

- normal stereo only

 fpm run --profile release -- -ns  file_name


## mixed block 
デフォルトでは long block, short block の他に、これらの二つの折衷で、低周波バンドは long block 高周波バンドは short block を利用する mixed block も使うようになっていますが、mixed block に対応していない再生機器やソフトがあります、この時 mixed block に切り替わる所でガザガザという音が少しします。そういう場合には mixed block を使わず long block と short block しか使わないようにするオプション -s を利用してください。

## 実行例

- debug run

 fpm run -- file_name

- release run

 fpm run --profile release -- file_name

- bitrate 320kBps, cut off frequency 32th band, normal stereo, no mixed blocks

 fpm run --profile release -- -b 14 -cut 32 -ns -s file_name