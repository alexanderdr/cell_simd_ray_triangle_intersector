/* This function casts a single ray against four triangles at once, taking advantage
   of the SPU's vector processing abilities.  It more than doubles the speed of this
   function, but there is some strange behavior when actually loading four triangles into
   a correctly aligned structure-- it appears to be free, because doing it in the function
   vs. pre-loading them once the triangles are in the local-store doesn't lead to any
   performance improvement.  I assume that it's free because of the dual-pipeline nature
   of the SPU and some strange behavior with the compiler.  The tricky part of converting
   this function from scalar to vector code is that when the ray misses the triangle in
   question, it can be known before the end of the function so the function can then
   return to skip unnecessary processing.  It uses a vector called "fail" to determine
   when one element of a vector has failed the check.  When all elements fail the function
   then returns.  Note that all vectors have to be 16 byte aligned or the program will fail
   during run time.
*/

/*
   The function of some of the SPU intrinsics used:
   spu_splats(float): returns a vector of floats filled with the argument
   spu_compeq(vector,vector): compares two vectors, if all bits are equal, the result is
                              all 1s, if not the result is all 0s.
   spu_sel(vector,vector,vector): returns a vector based on the bits of the
                                  third vector, if a bit is 0, then it has the same
                                  value as the first vector, if not, then the second
                                  vector.
   spu_gather(vector): the least significant bit of each element in the vector is gathered
                      then returned.  The result depends on the number of elements in the
                      vector.
   spu_extract(vector,int): Returns the element of the vector specified in the second
                            argument, used with gather to test for failure and to possibly
                            return early.  Int argument is usually 0 in this method.
   
*/
int cast4(Triangle* _tri, Ray* r,Tri_4 tri){
        
	Vec_4 p;
	Vec_4 s;
	Vec_4 q;
	Vec temp1;
	Vec temp2;
	vector float a, f, u, v, t;
	vector unsigned int fail;

	vector float zeroes = spu_splats(0.0f);

	vector unsigned int d;
	Vec pos;

	p = cross1_4(r->ray,tri.e2);
	a = dot4(tri.e1,p);
	vector unsigned int afix;
	fail = spu_cmpeq(a,zeroes);

	//__builtin_expect give the compiler a clue that it should generate prepare-to-branch assembly
	//for a particular choice which dramatically reduces the cost of branching if the prepared branch
	//is chosen.        
	if(__builtin_expect(spu_extract(spu_gather(fail),0)==15,0)){ //this adds about .4fps out of 10.4fps
		//puts("Joy");                         //even though the branch is never taken with "real" objects
		return -1;                           //The compiler is strange-- this probably causes dual-issueing
	}
	//if(a==0){ //scalar version of above statement
	//	return -1;
	//}
	
	
	//The following replaces the intermediate value of any triangles that are
	//perpendicular to the ray with something that won't cause division by
	//zero in the next step.  It is still known that they failed earlier
	//because of the fail vector.
	vector float ones = spu_splats(1.0f);
	vector float not_zeroes = {4321.0f,4321.0f,4321.0f,4321.0f};
	a = spu_sel(a,not_zeroes,fail);        
	

	f = recipf4(a);//1.0f/a;

	s = sub4(vec1to4(r->origin),tri.p1);//sub(&r->origin,&tri->p1,&s);
	u = spu_mul(f,dot4(s,p));

	//if(u < 0.0f || u > 1.0f) return -2;
	
	d = spu_or(spu_cmpgt(u,ones),spu_cmpgt(zeroes,u));
	fail = spu_or(fail,d);
	if(spu_extract(spu_gather(fail),0)==15) return -2; //15 == 1111, ie. four fails

	q = cross4(s,tri.e1);
	v = spu_mul(f,dot4(vec1to4(r->ray),q));

	d = spu_or(spu_cmpgt(zeroes,v),spu_cmpgt(spu_add(u,v),ones));
	fail = spu_or(fail,d);     
	if(spu_extract(spu_gather(fail),0)==15) return -3; //four fails
	//if(v < 0.0f || u + v > 1.0f) return -3;

	t = spu_mul(f,dot4(tri.e2,q));
	
	d = spu_cmpgt(zeroes,t);
	fail = spu_or(fail,d);

	int gather = spu_extract(spu_gather(fail),0);
	
	if(gather!=15){ //if there aren't 4 fails
	  
		//The follows calculates the number of failed valid triangles out of the 4
		unsigned int valid = 4-spu_extract(spu_cntb(((vector unsigned char)spu_splats(gather))),0);
		
		Triangle hit;
		unsigned int hitIndex;
		float hitT = 3487534.0f; //arbitrarily large number
		if(valid!=3){ //if valid = 3, then there were 3 failures, which means there's only 1 valid triangle
					  //if it's not, we have to iterate through the vector to find the one working triangle
			float tempT;
			unsigned int tempIndex = 0;
			while(tempIndex<4){
				if(spu_extract(fail,tempIndex)==-1){
					//printf("Skip\n");
					++tempIndex;
					continue;
				}
				tempT = spu_extract(t,tempIndex);
				if(tempT<hitT){
					//printf("%x\n",spu_extract(fail,tempIndex));
					hitT = tempT;
					hitIndex = tempIndex;
				}
				++tempIndex;
			}
			hit = _tri[hitIndex];
			//printf("Complex ");
		} else {
			hitIndex = 31-spu_extract(spu_cntlz(((vector unsigned int)spu_splats(gather))),0);
			hit = _tri[hitIndex];//extractTriangle(tri,hitIndex);
			hitT = spu_extract(t,hitIndex);
			//printf("Simple ");
		}


		if(spu_extract(fail,hitIndex)==-1) printf("Said there's a hit but no hit\n");
		
		if(fabs(hitT)>fabs(r->t)) return 2;
		
		temp1 = scaleNewVec(&hit.e1,spu_extract(u,hitIndex));
		temp2 = scaleNewVec(&hit.e2,spu_extract(v,hitIndex));
		add(&temp1,&temp2,&pos);
		addInPlace(&pos,&hit.p1);
		r->intersection.x = pos.x;
		r->intersection.y = pos.y;
		r->intersection.z = pos.z;
		r->itri = &_tri[hitIndex];//&hit;
		r->color.x = .3;
		r->color.y = .7;
		r->color.z = 1;
		r->t = hitT;
		return 1;
	}

	return -4;

}