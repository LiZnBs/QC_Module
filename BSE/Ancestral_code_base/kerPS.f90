pairkk = k**2;pairkq = k*q*(zk*zq + y*Sqrt(1 - zk**2)*Sqrt(1 - zq**2));pairkP = (0,1)*k*P*zk;pairqq = q**2;pairPq = (0,1)*P*q*zq;pairPP = -P**2;select case(nKerPS*i+j)    case(1)        Kernel =         (16*GLPP*(4*Bms*Bpu - Ams*Apu*pairPP + 4*Ams*Apu*pairqq))/&
          ((Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(2)        Kernel =         (32*GLPP*((Apu*Bms + Ams*Bpu)*pairPP + 2*(Apu*Bms - Ams*Bpu)*pairPq))/&
          ((Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(3)        Kernel =         (32*GLPP*pairPq*((Apu*Bms + Ams*Bpu)*pairPq + 2*(Apu*Bms - Ams*Bpu)*pairqq))/&
          ((Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(4)        Kernel =         (64*Ams*Apu*GLPP*(pairPq**2 - pairPP*pairqq))/&
          ((Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(5)        Kernel =         (32*GLPP*(-(pairkk**2*((Apu*Bms + Ams*Bpu)*pairPP + 2*(Apu*Bms - Ams*Bpu)*pairPq)) + &
              pairkk*((Apu*Bms + Ams*Bpu)*pairkP**2 - 2*Apu*Bms*pairPq**2 - 2*Ams*Bpu*pairPq**2 + &
                 2*pairkq*((Apu*Bms + Ams*Bpu)*pairPP + 4*(Apu*Bms - Ams*Bpu)*pairPq) + &
                 2*pairkP*((Apu*Bms - Ams*Bpu)*pairkq + (Apu*Bms + Ams*Bpu)*pairPq) - &
                 Apu*Bms*pairPP*pairqq - Ams*Bpu*pairPP*pairqq - 6*Apu*Bms*pairPq*pairqq + &
                 6*Ams*Bpu*pairPq*pairqq) + &
              pairkP*(-((Apu*Bms + Ams*Bpu)*pairkP*(4*pairkq - pairqq)) + &
                 2*pairkq*((-4*Apu*Bms + 4*Ams*Bpu)*pairkq + (Apu*Bms + Ams*Bpu)*pairPq + &
                    3*(Apu*Bms - Ams*Bpu)*pairqq))))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(6)        Kernel =         (16*GLPP*(pairkk**2*(-(Ams*Apu*pairPP**2) + 8*Ams*Apu*pairPq**2 + &
                 4*pairPP*(Bms*Bpu - Ams*Apu*pairqq)) + &
              pairkP*(2*pairkq*pairPq*(-4*Bms*Bpu + 16*Ams*Apu*pairkq + Ams*Apu*pairPP - &
                    8*Ams*Apu*pairqq) - pairkP*(4*pairkq - pairqq)*&
                  (-4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq)) + &
              pairkk*(8*Bms*Bpu*pairPq**2 - 2*Ams*Apu*pairPP*pairPq**2 + 4*Bms*Bpu*pairPP*pairqq - &
                 Ams*Apu*pairPP**2*pairqq + 16*Ams*Apu*pairPq**2*pairqq - &
                 4*Ams*Apu*pairPP*pairqq**2 + &
                 pairkP**2*(-4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) + &
                 2*pairkP*pairPq*(-4*Bms*Bpu - 4*Ams*Apu*pairkq + Ams*Apu*pairPP + &
                    4*Ams*Apu*pairqq) + pairkq*&
                  (2*Ams*Apu*pairPP**2 - 32*Ams*Apu*pairPq**2 + &
                    pairPP*(-8*Bms*Bpu + 8*Ams*Apu*pairqq)))))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(7)        Kernel =         (16*GLPP*pairPq*(pairkk**2*pairPq*(4*Bms*Bpu - Ams*Apu*pairPP + 4*Ams*Apu*pairqq) + &
              pairkP*(2*Ams*Apu*pairkP*pairPq*pairqq + &
                 4*pairkq**2*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) - &
                 pairkq*(8*Ams*Apu*pairkP*pairPq - 4*Ams*Apu*pairPq**2 + &
                    3*pairqq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq))) + &
              pairkk*(2*Ams*Apu*pairkP**2*pairPq - &
                 pairkP*(-4*Ams*Apu*pairPq**2 + &
                    pairkq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq)) + &
                 pairPq*(-4*Ams*Apu*pairPq**2 - 16*pairkq*(Bms*Bpu + Ams*Apu*pairqq) + &
                    pairqq*(12*Bms*Bpu + Ams*Apu*pairPP + 12*Ams*Apu*pairqq)))))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(8)        Kernel =         (32*GLPP*(2*(Apu*Bms + Ams*Bpu)*pairkk**2*(pairPq**2 - pairPP*pairqq) + &
              pairkP*(4*pairkq**2*((Apu*Bms - Ams*Bpu)*pairPP + 2*(Apu*Bms + Ams*Bpu)*pairPq) + &
                 pairkP*pairqq*((Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq) - &
                 pairkq*((-2*Apu*Bms + 2*Ams*Bpu)*pairPq**2 + 3*(Apu*Bms - Ams*Bpu)*pairPP*pairqq + &
                    2*(Apu*Bms + Ams*Bpu)*pairPq*pairqq + &
                    4*pairkP*((Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq))) + &
              pairkk*(pairkP**2*((Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq) - &
                 2*((Apu*Bms - Ams*Bpu)*pairPq - (Apu*Bms + Ams*Bpu)*pairqq)*&
                  (pairPq**2 - pairPP*pairqq) + &
                 pairkP*(-(pairkq*((Apu*Bms - Ams*Bpu)*pairPP + 2*(Apu*Bms + Ams*Bpu)*pairPq)) + &
                    2*pairPq*((Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq)) + &
                 2*pairkq*(-4*(Apu*Bms + Ams*Bpu)*pairPq**2 + &
                    pairPP*((-(Apu*Bms) + Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq)))))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(9)        Kernel =         (64*GLPP*((Apu*Bms + Ams*Bpu)*pairkP**3 + &
              2*pairkP**2*((Apu*Bms - Ams*Bpu)*pairkq - (Apu*Bms + Ams*Bpu)*pairPq + &
                 (-(Apu*Bms) + Ams*Bpu)*pairqq) + &
              pairkP*(-(pairkk*((Apu*Bms + Ams*Bpu)*pairPP + (-(Apu*Bms) + Ams*Bpu)*pairPq)) + &
                 pairkq*((Apu*Bms + Ams*Bpu)*pairPP + 4*(-(Apu*Bms) + Ams*Bpu)*pairPq) + &
                 pairPq*((Apu*Bms + Ams*Bpu)*pairPq + 3*(Apu*Bms - Ams*Bpu)*pairqq)) + &
              pairPP*(pairkk*((-3*Apu*Bms + 3*Ams*Bpu)*pairkq + (Apu*Bms + Ams*Bpu)*pairPq + &
                    2*(Apu*Bms - Ams*Bpu)*pairqq) + &
                 pairkq*(4*(Apu*Bms - Ams*Bpu)*pairkq - (Apu*Bms + Ams*Bpu)*pairPq + &
                    3*(-(Apu*Bms) + Ams*Bpu)*pairqq))))/&
          (3.*pairkP*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(10)        Kernel =         (32*GLPP*(2*pairkP**2*(4*Bms*Bpu - 4*Ams*Apu*pairkq - Ams*Apu*pairPP)*pairPq + &
              pairkP**3*(-4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) + &
              pairPP*pairPq*(pairkk*(-4*Bms*Bpu + 12*Ams*Apu*pairkq + Ams*Apu*pairPP - &
                    4*Ams*Apu*pairqq) + pairkq*&
                  (4*Bms*Bpu - 16*Ams*Apu*pairkq - Ams*Apu*pairPP + 8*Ams*Apu*pairqq)) + &
              pairkP*(pairPq**2*(-4*Bms*Bpu + Ams*Apu*pairPP - 8*Ams*Apu*pairqq) + &
                 pairkk*(-(Ams*Apu*pairPP**2) - 4*Ams*Apu*pairPq**2 + &
                    4*pairPP*(Bms*Bpu - Ams*Apu*pairqq)) + &
                 pairkq*(Ams*Apu*pairPP**2 + 16*Ams*Apu*pairPq**2 + &
                    pairPP*(-4*Bms*Bpu + 4*Ams*Apu*pairqq)))))/&
          (3.*pairkP*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(11)        Kernel =         (16*GLPP*pairPq*(4*Ams*Apu*pairkP**3*pairPq - &
              2*pairkP**2*(4*Ams*Apu*pairPq**2 + &
                 pairkq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) - &
                 pairqq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq)) - &
              pairkP*pairPq*(-4*Ams*Apu*pairPq**2 + 12*Bms*Bpu*pairqq + 3*Ams*Apu*pairPP*pairqq + &
                 12*Ams*Apu*pairqq**2 - 8*pairkq*(2*Bms*Bpu + Ams*Apu*pairPP + 2*Ams*Apu*pairqq) + &
                 pairkk*(4*Bms*Bpu + 5*Ams*Apu*pairPP + 4*Ams*Apu*pairqq)) + &
              pairPP*(pairkk*(4*Ams*Apu*pairPq**2 + &
                    3*pairkq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) - &
                    2*pairqq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq)) + &
                 pairkq*(-4*Ams*Apu*pairPq**2 - &
                    4*pairkq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq) + &
                    3*pairqq*(4*Bms*Bpu + Ams*Apu*pairPP + 4*Ams*Apu*pairqq)))))/&
          (3.*pairkP*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(12)        Kernel =         (32*GLPP*(2*pairkP**3*((Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq) - &
              2*pairkP**2*(2*(Apu*Bms - Ams*Bpu)*pairPq**2 + &
                 pairkq*((Apu*Bms - Ams*Bpu)*pairPP + 2*(Apu*Bms + Ams*Bpu)*pairPq) + &
                 (-(Apu*Bms) + Ams*Bpu)*pairPP*pairqq + 2*(Apu*Bms + Ams*Bpu)*pairPq*pairqq) - &
              pairkP*(pairPq*((-2*Apu*Bms + 2*Ams*Bpu)*pairPq**2 + &
                    3*(Apu*Bms - Ams*Bpu)*pairPP*pairqq + 2*(Apu*Bms + Ams*Bpu)*pairPq*pairqq) - &
                 2*pairkq*(4*(Apu*Bms + Ams*Bpu)*pairPq**2 + &
                    pairPP*(3*(Apu*Bms - Ams*Bpu)*pairPq + 2*(Apu*Bms + Ams*Bpu)*pairqq)) + &
                 pairkk*(2*(Apu*Bms + Ams*Bpu)*pairPq**2 + &
                    pairPP*(3*(Apu*Bms - Ams*Bpu)*pairPq + 4*(Apu*Bms + Ams*Bpu)*pairqq))) + &
              pairPP*(pairkq*((-2*Apu*Bms + 2*Ams*Bpu)*pairPq**2 - &
                    4*pairkq*((Apu*Bms - Ams*Bpu)*pairPP + 2*(Apu*Bms + Ams*Bpu)*pairPq) + &
                    3*(Apu*Bms - Ams*Bpu)*pairPP*pairqq + 2*(Apu*Bms + Ams*Bpu)*pairPq*pairqq) + &
                 pairkk*(3*pairkq*((Apu*Bms - Ams*Bpu)*pairPP + 2*(Apu*Bms + Ams*Bpu)*pairPq) + &
                    2*(Apu*Bms - Ams*Bpu)*(pairPq**2 - pairPP*pairqq)))))/&
          (3.*pairkP*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(13)        Kernel =         (64*Ams*Apu*GLPP*(2*pairkP**2*pairqq + pairkq*pairPP*pairqq - &
              pairkP*pairPq*(2*pairkq + pairqq) + &
              pairkk*(pairkq*pairPP - pairkP*pairPq + 2*pairPq**2 - 2*pairPP*pairqq)))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(14)        Kernel =         (-64*(Apu*Bms + Ams*Bpu)*GLPP*(2*pairkP**2*pairqq + pairkq*pairPP*pairqq - &
              pairkP*pairPq*(2*pairkq + pairqq) + &
              pairkk*(pairkq*pairPP - pairkP*pairPq + 2*pairPq**2 - 2*pairPP*pairqq)))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(15)        Kernel =         (32*(Apu*Bms - Ams*Bpu)*GLPP*pairPq*&
            (2*pairkP**2*pairqq + pairkq*pairPP*pairqq - pairkP*pairPq*(2*pairkq + pairqq) + &
              pairkk*(pairkq*pairPP - pairkP*pairPq + 2*pairPq**2 - 2*pairPP*pairqq)))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))    case(16)        Kernel =         (16*GLPP*(4*Bms*Bpu + Ams*Apu*pairPP - 4*Ams*Apu*pairqq)*&
            (-2*pairkP**2*pairqq - pairkq*pairPP*pairqq + pairkP*pairPq*(2*pairkq + pairqq) + &
              pairkk*(-(pairkq*pairPP) + pairkP*pairPq - 2*pairPq**2 + 2*pairPP*pairqq)))/&
          (3.*(pairkP**2 - pairkk*pairPP)*(pairkk - 2*pairkq + pairqq)*&
            (Ams**2*pairPP + 4*(Bms**2 - Ams**2*pairPq + Ams**2*pairqq))*&
            (Apu**2*pairPP + 4*(Bpu**2 + Apu**2*pairPq + Apu**2*pairqq)))end select