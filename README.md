# uzura3_fpm

MP3 encoder in Fortran95.

command line 引数を取る所と stream access 以外は（これらは Fortran 2003 になります)、Fortran95 で書いてあると思います。

15 年放置の後、形式的エンコーダー部分のバグが取れました。その結果、心理音響部のバグも取ることが出来ました。

思いのほか良い音が出るようになりました。

## bug 

mp3 のプログラムは scale factor を探すところが高次ベクトル空間の残差最小化の問題と見なせます。この部分に計算に時間がかかっていたので早めに探索を諦める if 文を１個入れていたのですが、それの筋が悪くて悪さをしていました。最初 dual Cerelon 300A 機で作り始めて実行がとても遅かったので、思い付きで入れたのですがすっかり忘れていました･･･

これを外した後は、理屈通りにパラメータ変化に従って動くようになり、音を聞くことで心理音響解析部の誤りも修正できました。

心理音響部については、そもそも dB とか音響学の基礎レベルの理解が全くできておらず音圧の power と absolute amplitude を混同して物理量の次元が狂っているなど、プログラム以前に内容がひどく間違っていました。ただ基本的なプログラム構造は大丈夫でサブルーチン内の修正で済んだのが不幸中の幸いでした。

## Psychoacoustic model  

心理音響解析のモデルは Bosse Lincoln の "An Experimental High Fidelity Perceptual Audio Coder Project in MUS420 Win 97" (1998) におおむね従っています。ただマスキングの効果は線形和ではなく非線形和を取っています。これは uzura mp1/mp2 encoder と同様です。 論文中でも非線形和について触れられているのですが、モデルとしては採用されていません。なお Lincoln 論文中でもマスク和のところで power と absolute amplitude が混乱しているように見えるので、適宜修正してプログラムしました。

非線形和（ノルムの定義）に関しては、F. Baumgarte, C. Ferekidis and H. Fuchs の "A Nonlinear Psychoacoustic Model Applied to the ISO MPEG Layer 3 Coder" も参考にしました。非線形和の指数は、こちらに従って 0.3 にしてあります。これと組み合わせる Spreading function も Lincolnn のものではなく、単純な線形のものを取っています。

## bit distribution

心理音響解析によりマスキングの閾値が求められた後、mp3 フォーマット では 1 frame 中の granule/channel 間の bit 分配がとても重要になります。上記 Lincoln の論文は独自のフォーマットを使っているので、この部分は uzura 独自のものです。かなり適当に作りましたが音質にとって極めて重要な部分で、ここが悪いと scale factor を回しても音質の向上余地が限られます。

## mixed block 

デフォルトでは long block, short block の他に、これらの二つの折衷で低周波バンドは long block、 高周波バンドは short block を利用する mixed block も使うようになっています。この点が Uzura の特徴と言えるでしょうか。しかしこれは音質向上にさほど寄与しないようです。

ところで再生機器やソフトによっては、 mixed block に対応していないものがあります。この場合 mixed block に切り替わる所でガサガサゴボゴボという音が少しします。そういう時には mixed block を使わないようにするオプション -s を利用してください。

SONY Walkman NW A17 は対応、SONY Xperia 8 は非対応でした。

## Mid Side Stereo

Normal Stereo と Mid-Side stereo の切り替えアルゴリズムが素朴すぎて、うまくいかない時があります。（昔のアナログ時代によくあった左右のチャンネルに別々の楽器を Mix したステレオなど。）この場合は、下記の実行例を参考に Normal Stereo にするか、Normal Stereo になりやすいパラメータを選んで下さい。

## 実行例 fortran package manager (fpm) 利用例

実行は Fortran Package Manager (fpm) を利用すると簡単にできます。　https://github.com/fortran-lang/fpm

心理音響解析のマスキングの畳み込みの非線形和にとても時間がかかるので、Haswell i7 CPU 3.4GHz だと、エンコードするのに曲 1 分当たり 4 分くらいかかります。心理音響解析部は計算精度もいらないし演算ヘビーでランダムなメモリーアクセスもないので高速化の余地は大いにあると思います。

** file_name.wa --> file_name.mp3 **

- debug run

 fpm run -- file_name

- release run

 fpm run --profile release -- file_name

- run with options: bitrate 320kBps, normal stereo, no mix-block s

 fpm run --profile release -- -b 14 -ns -s file_name

- non-linear sum; p-norm for adding spreading function  

 fpm run --profile release -- -pow 0.5 -- file_name 

- lineaar sum; 1-norm

 fpm run --profile release -- -pow 1.0 -- file_name 

- lineaar sum; 1-norm (another formulation)

 fpm run --profile release -- -1norm -- file_name 

 - enhance tonality

 fpm run --profile release -- -tonality 2.0 -- file_name 

### Mid-Side stereo

 normal stereo と mid-side の切り替えがうまく行かない時があります。オプション xms 小で Normal Stereo、xms 大で Mid-Side になりやすくなります

- normal stereo 多め

 fpm run --profile release -- -xms 0.2  file_name

- normal stereo only

 fpm run --profile release -- -ns  file_name
