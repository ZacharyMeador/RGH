//This is just a block from the entire code



for(int i=0; i<ndchit; i++) {
	    int it=(*hitn)[i]-1;
        
//        cout << it << " x-component = " << (*vx)[it]<< " and " << " y-component = " << (*vy)[it] << endl;
         const double pi = 3.14159265358979323846;
         int nsect;
         double phi;
         int x = (*vx)[it];
         int y = (*vy)[it];
         phi = atan2(y,x)*180/pi;
//         cout << it << " has an angle of "<<phi << endl;
            if(0<y){
//                  cout << "y-component is positive " << y << endl;
                 if(0>x){
//                      cout << "x-component is negative " << (*vx)[it] << endl;
                     phi = phi + 180;
                 // cout << it << "'s angle in degrees = " << phi << endl;
                 }
                else {
                    phi = phi;
                }
              }
             else if(0>y){
                 if(0>x){
                 phi = phi + 180;
//                  cout << it << "'s angle in degrees = " << phi << endl;
                 }
                 else{
                 phi = phi + 360;
//                cout << it << "'s angle in degrees = " << phi << endl;
                 }
             }
            if(30<phi && phi<=90){
                nsect = 2;
                count_sect2 += 1;
//                 cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
            }
            else if(90<phi && phi<=150){
                nsect = 3;
                count_sect3 += 1;
//                 cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
            }
            else if(150<phi && phi<=210){
            nsect = 4;
            count_sect4 += 1;
//             cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
            }
            else if(210<phi && phi<=270){
            nsect = 5;
            count_sect5 += 1;
//             cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
            }
            else if(270<phi && phi<=330){
                nsect = 6;
                count_sect6 += 1;
//                 cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
            }
            else {
                nsect = 1;
                count_sect1 += 1;
//                 cout << it << " is in sector " << nsect << " with " << phi << " deg " << " x-component = " << (*vx)[it] << " y-component = "<<(*vy)[it]<< endl;
           }
           
   
