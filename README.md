# FUN-ECM2017K  
keifujunagaが作成したプログラム  
9/5  
scalar.cのコアダンプ（浮動小数点演算例外）の修正  
原因：scakar関数の引数を１つ多く設定していた^q^  

9/16  
●normal_add.c double_add.cのコメントの修正  
●point.cの関数が正しくない  
・protomon  
mongomery曲線のy座標を復元後、stage1のkPのy座標を復元  
・他2つ  
1-Yや-d-1など  

9/18  
●extend edwards曲線とmontgomery曲線の両方とも使えるようにしました  
●point.cのprotomon以外の修正  

9/27  
●point.cの修正  
●scalar.cのmscalar修正  
●normal_add.cのmontgomery_add修正  

10月  
●モントゴメリーのスカラー倍と座標同士の足し算のバグ修正

11/1  
●プログラムを用いて計算を開始  
●プログラムはモントゴメリーを使用して素数テーブルを実装したものを使用  


