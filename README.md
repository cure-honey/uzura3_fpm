# uzura3_fpm

MP3 encoder in Fortran95.

command line 引数を取る所と stream access 以外は（これらは Fortran 2003 になります)、Fortran95 で書いてあると思います。

15 年放置の後、scale factor loop のバグ（仕様？）が取れました。心理音響部は新たに作り直しました。

今のところパラメータ的に高周波の方にビットを割り振り過ぎて、少しキンキンする疲れる音になっています。

## bug 

mp3 のプログラムは scale factor を探すところが高次ベクトル空間の残差最小化の問題と見なせますが、ここに問題がありました。


## Psychoacoustic model  

心理音響解析のモデルは F. Baumgarte, C. Ferekidis and H. Fuchs の "A Nonlinear Psychoacoustic Model Applied to the ISO MPEG Layer 3 Coder" におおむね従っています。

## bit distribution

心理音響解析によりマスキングの閾値が求められた後、mp3 フォーマット では 1 frame 中の granule/channel 間の bit 分配がとても重要になります。scale factor を回すより大事かもしれません。

## mixed block 

デフォルトでは long block, short block の他に、mixed block も使うようになっています。この点が Uzura の特徴と言えるでしょうか。しかしこれは音質向上にさほど寄与しないようです。Mixed block は、低周波バンドは long block 高周波バンドは short block を利用する折衷になっています。

ところで再生機器やソフトによっては、 mixed block に対応していないものがあります。この場合 mixed block に切り替わる所でガサガサゴボゴボという音が少しします。そういう時には mixed block を使わないようにするオプション -s を利用してください。



## Mid Side Stereo

Normal Stereo と Mid-Side stereo の切り替えアルゴリズムが素朴すぎて、うまくいかない時がママあります。（昔のアナログ時代によくあった左右のチャンネルに別々の楽器を Mix したステレオなど。）下記の実行例を参考に Normal Stereo にするか、Normal Stereo になりやすいパラメータを選んで下さい。

## 実行例 fortran package manager (fpm) 利用例

実行は Fortran Package Manager (fpm) を利用すると簡単にできます。 https://github.com/fortran-lang/fpm

** file_name.wa --> file_name.mp3 **

- debug run

 fpm run -- file_name

- release run

 fpm run --profile release -- file_name

- run with options: bitrate 320kBps, normal stereo, no mix-block s

 fpm run --profile release -- -b 14 -ns -s file_name

- no scalefactor loop

 fpm run --profile release -- -inner file_name

### Mid-Side stereo

 normal stereo と mid-side の切り替えがうまく行かない時があります。オプション xms 小で Normal Stereo、xms 大で Mid-Side になりやすくなります

- normal stereo 多め

 fpm run --profile release -- -xms 0.2  file_name

- normal stereo only

 fpm run --profile release -- -ns  file_name
