import React, { useEffect, useState } from 'react';
import './swiper.css'
import Butterfly from '../../image/8.png'
import Butterfly9 from '../../image/9.png'
import Man1 from '../../image/1.png'
import Woman1 from '../../image/2.png'
import Man2 from '../../image/3.png'
import Woman2 from '../../image/4.png'
import Man3 from '../../image/6.png'
import Woman3 from '../../image/5.png'
import Eropen700 from '../../image/EROPEN700.png'
import { Swiper, SwiperSlide } from "swiper/react";
import { Parallax, Pagination, Navigation, Autoplay, EffectFade, Lazy } from "swiper";
import "swiper/css";
import "swiper/css/effect-fade";
import "swiper/css/navigation";
import "swiper/css/pagination";
import "swiper/css/parallax";
import "swiper/css/lazy";
import { useTranslation } from 'react-i18next';


const Swipers = () => {
    const { t } = useTranslation()
    const [slider] = useState([
        { id: 1, imageMan: Man1, imageWoman: Woman1, title: t('swiper.swiperText1')},
        { id: 2, imageMan: Man2, imageWoman: Woman2, title: t('swiper.swiperText2')},
        { id: 3, imageMan: Woman3, imageWoman: Man3, title: t('swiper.swiperText3')}
    ])
    const [ active, setActive] = useState(false)
    useEffect(() => {
        setActive(true)
    },[])

    return (
        <div className={active? "main-content --anim-items --animate" : 'main-content --anim-items'}>
            <div className="main-content__row">
                <div className='main-content__slider'>
                    <Swiper
                        modules={[EffectFade, Parallax, Pagination, Navigation, Autoplay, Lazy]}
                        navigation={true}
                        spaceBetween={30}
                        speed={1500}
                        lazy={true}
                        // parallax={true}
                        loop={true}
                        slideToClickedSlide={true}
                        effect={"fade"}
                        fadeEffect = {{
                             crossFade: true }}
                        grabCursor={true}
                        keyboard = {{
                            enabled: true,
			                onlyInViewport: true,
			                pageUpDown: true
                        }}
                        // preloadImages =  {!1},
                        autoplay={{
                            delay: 2500,
                            disableOnInteraction: false,
                        }}
                        pagination={{
                            el: '.swiper-pagination',
                            clickable: true,
                            dynamicBullets: true
                        }}
                    >
                        {slider.map(slide => (
                            <SwiperSlide key={slide.id}>
                                {({ isActive }) => (
                                    <div className="main-content__image">
                                        <div className="main-content__brand main-content__brand">
                                            <span className="main-content__subtitle">
                                                <div className="brand__logo" style={
                                                    isActive ? {
                                                        
                                                        transitionDuration: '2000ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(0px, 0px, 0px) scale(1)'
                                                        
                                                        
                                                    } : 
                                                    {
                                                        transitionDuration: '2000ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(0px, 0px, 200px) scale(3)'
                                                    }
                                                }>
                                                    <img src={Eropen700} alt="EROPEN-700" />
                                                </div>
                                            </span>
                                        </div>
                                        <div className="main-content__people --anim-items --animate">
                                            <div className="main-content__people-image"
                                                style={
                                                    isActive ? {
                                                        transitionDuration: '1500ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(0px, 0px, 0px)'
                                                    } : {
                                                        transitionDuration: '1500ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(300px, 0px, 0px)'
                                                    }
                                                }
                                            >
                                                <img src={slide.imageMan} alt="Man" />
                                                <div className="main-content__people-flym">
                                                    <img src={Butterfly} alt="Butterfly" />
                                                </div>
                                            </div>
                                            <div className="main-content__people-image"
                                                style={
                                                    isActive ? {
                                                        transitionDuration: '1500ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(0px, 0px, 0px)'
                                                    } : {
                                                        transitionDuration: '1500ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(-300px, 0px, 0px)'
                                                    }
                                                }
                                            >
                                                <img src={slide.imageWoman} alt="Woman" />
                                                <div className="main-content__people-flyw">
                                                    <img src={Butterfly9} alt="Butterfly" />
                                                </div>
                                            </div>
                                        </div>
                                        <h2 className="main-content__title" style={
                                                    isActive ? {
                                                        transitionDuration: '2000ms',
                                                        opacity: 1,
                                                        transform: 'translate3d(0px, 0px, 0px) scale(1)'
                                                        
                                                    } : {
                                                        transitionDuration: '0ms',
                                                        opacity: 0,
                                                        transform: 'translate3d(px, 0px, 300px) scale(0)'
                                                    }
                                                }>
                                            {slide.title}
                                        </h2>
                                    </div>
                                )}
                            </SwiperSlide>
                        ))}
                    </Swiper>
                </div>
            </div>
        </div>
    );
};

export default Swipers;