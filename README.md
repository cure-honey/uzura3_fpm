# uzura3_fpm

MP3 encoder in Fortran.

すこしづつ Modern Fortran に直したいです。
とりあえず小文字にして fpm + gfortran でコンパイル出来るようにしました。 

今の所中身は２０年前のママです。
心理音響解析のモデルが良くないので、mp1/mp2 に採用したものと同じ原理に直したいと思っているのですが･･･
（周波数と強度を非線形変換した後、線形演算する。n-norm; n 非整数, n<1）

## mixed block 
デフォルトでは long block, short block の他に、これらの二つの折衷で、低周波バンドは long block 高周波バンドは short block を利用する mixed block も使うようになっていますが、mixed block に対応していない再生機器やソフトがあります、この時 mixed block に切り替わる所でガザガザという音が少しします。そういう場合には mixed block を使わず long block と short block しか使わないようにするオプション -s を利用してください。

## 実行例

- debug run

 fpm run -- file_name

- release run

 fpm run --profile release -- file_name

- bitrate 320kBps, cut off frequency 32 band, normal stereo, no mix-block s

 fpm run --profile release -- -b 14 -cut 32 -ns -s file_name