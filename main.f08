program Excer_5
    implicit none
    real(8), dimension(:), allocatable:: dens, fi, E
    real(8), dimension(:,:), allocatable:: x_ion, p_ion, x_e, p_e
    real(8):: L, d, T, dt, dx, q_e, q_p, m_e, m_p, dL, E_0, B_max, B, c, R_m, B_min, temp
    integer(8):: ne_0, N_g, i, j, n_c, quantinty,  quantinty_change, counter_last

    ne_0 = 8E10                            !Концентрация электронов
    L = 200                                !Длина просраносва
    d = 1.2                                !Ширина слоя
    N_g = 20000                            !Количесво ячейк
    T = 5E-8                               !Время моделирования
    dt = T/N_g                             !Шаг по времени
    dx = L/N_g                             !Шаг по расстоянию
    q_e = -4.8E-10                         !Заряд электрона
    q_p = 4.8E-10                          !Заряд протона
    m_e = 9.1E-28                          !масса электрона
    m_p = 1.7E-24                          !масса протона
    n_c = 1000                             !Количесво счетных
    E_0 = 1.8E-5                           !Нач. энергия(Эрг)
    B_max = 5000                           !Индукция (Гс)
    B_min = 4500
    c = 3E10
    R_m = m_p/m_e
    temp = (d/dx)/2


    allocate(dens(0:N_g + 1))
    allocate(fi(0:N_g + 1))
    allocate(E(0:N_g + 1))
    allocate(x_ion(0:n_c,0:N_g + 1))           !коор йонов. 1-номер йона, 2-момент времени
    allocate(p_ion(0:n_c,0:N_g + 1))           !импульс йонов влодь Х. 1-номер йона, 2-момент времени
    allocate(x_e(0:n_c,0:N_g + 1))             !коор е. 1-номер е, 2-момент времени
    allocate(p_e(0:n_c,0:N_g + 1))             !импульс е влодь Х. 1-номер е, 2-момент времени


    counter_last = 0
    quantinty = n_c

    fi(0) = 0
    fi(1) = 0
    fi(N_g) = 0
    E(0) = 0
    E(N_g) = 0


    do i = 0, n_c                                                                   !Нач. скорость элек и ионов
        p_ion(i,0) = 0
        p_e(i,0) = 0
    enddo

    dL = (-1) * d/2
    i = 0
    do while(quantinty > 0)
        if(i >= int((N_g/2)-60) .and. i <= int((N_g/2)+60) .and. quantinty > 0) then    !Расп. электронов и йонов           
            quantinty_change = int(quantinty * abs(cos(3.1415 * abs(dL)/d)))
            if(quantinty_change == 0) then
                quantinty_change = 1
            endif
            
            
            if(quantinty <= quantinty_change) then

                do j = counter_last, counter_last + quantinty
                    x_e(j,0) = i * dx
                    x_ion(j,0) = i * dx
                enddo

                exit

            else
                do j = counter_last, counter_last + quantinty_change
                    x_e(j,0) = i * dx
                    x_ion(j,0) = i * dx
                enddo
            endif

            counter_last = counter_last + quantinty_change
            quantinty = quantinty -  quantinty_change
            dL = dL + dx
        endif

        E(i) = 0
        fi(i) = 0
        dens(i) = 0

        if(i >= N_g - 1) then
            i = 0
            dL = (-1) * d/2
        else
            i = i + 1
        endif

    enddo


    do j = 0, N_g - 1                                                                       !Цикль по времени

        do i = 0, n_c

            dens(int(x_ion(i,j)/dx)) = dens(int(x_ion(i,j)/dx)) + ((ne_0/n_c) * q_p)
            dens(int(x_e(i,j)/dx)) = dens(int(x_e(i,j)/dx)) + ((ne_0/n_c) * q_e) 

        enddo


        do i = 0, N_g
            fi(0) = fi(0) + i*dens(i)

        enddo

        fi(0) = fi(0)/N_g
        fi(1) = dens(0) + 2*fi(0)
        

        do i = 2, N_g
            fi(i) = (-1)*dens(i - 1)*dx**2 + 2*fi(i - 1) - fi(i - 2)
        enddo

        do i = 1, N_g - 2
            E(i) = (fi(i - 1) - fi(i + 1))/(2 * dx)
        enddo

        do i = 0, n_c                                                                       !Цикль по йонам и электроном
            
            if(x_e(i,j) >= L) then
                x_e(i,j) = L - dx
            else if(x_e(i,j) <= 0) then
                x_e(i,j) =  dx
            endif

            if(x_ion(i,j + 1) >= L) then
                x_ion(i,j) = L - dx
            else if(x_ion(i,j) <= 0) then
                x_ion(i,j) =  dx
            endif



            B = ((B_min - B_max)/(100.0))  * x_e(i,j) + ((B_min - B_max)/200.0) * L + B_max
            p_e(i,j + 1) = p_e(i,j) + (((-1) * (E_0/B) * ((B_min - B_max)/100.0)) + (q_e * E(int(x_e(i,j)/dx)))) * dt
            if(p_e(i,j) >= 0) then
                x_e(i,j + 1) = x_e(i,j) + sqrt(1/(1/(c**2) + (m_e**2)/(p_e(i,j)**2))) * dt
            else
                x_e(i,j + 1) = x_e(i,j) - sqrt(1/(1/(c**2) + (m_e**2)/(p_e(i,j)**2))) * dt
            endif



            p_ion(i,j + 1) = p_ion(i,j) + (q_p * (E(int(x_ion(i,j)/dx))/R_m)) * dt
            if(p_ion(i,j) >= 0) then
                x_ion(i,j + 1) = x_ion(i,j) + sqrt(1/(1/(c**2) + (m_p**2)/(p_ion(i,j)**2))) * dt
            else
                x_ion(i,j + 1) = x_ion(i,j) - sqrt(1/(1/(c**2) + (m_p**2)/(p_ion(i,j)**2))) * dt
            endif




        enddo


        do i = 0, N_g
            dens(i) = 0
        enddo

    enddo

    open(1, file = "spectr_E.csv")
    open(2, file = "spectr_Ion.csv")
    open(3, file = "coor.csv")
    open(4, file = "last_I.csv")

    do i = 0, n_c
        write(1,*) i,",",x_e(i,N_g-1),",", ((abs(p_e(i,N_g-1))*(c**2))/(sqrt(1/(1/(c**2)&
         + (m_e**2)/(p_e(i,N_g-1)**2))))) *6.3E12
        write(2,*) i,",",x_ion(i,N_g-1),",",((abs(p_ion(i,N_g-1))*(c**2))/(sqrt(1/(1/(c**2)&
         + (m_p**2)/(p_ion(i,N_g-1)**2)))))  *6.3E12
         write(4,*) i, ",", x_ion(i,0), ",", x_ion(i,18400),",",x_ion(i,N_g - 1)
    enddo

    do i = 0, N_g
        write(3,*) i * dt,",",x_ion(0,i),",",x_ion(int(n_c/2),i),",",x_ion(int(n_c-1),i),&
        ",",x_e(0,i),",",x_e(int(n_c/2),i),",",x_e(n_c-1,i)
    enddo


    close(1)
    close(2)
    close(3)
    print*, "Done!"
end program Excer_5