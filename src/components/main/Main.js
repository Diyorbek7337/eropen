import React, { useEffect, useState } from 'react';
import './main.css'
// import Logo from '../../image/logo.png';
import "@google/model-viewer/dist/model-viewer";
import Popup from '../../image/popup-image.png'
import Eropen700 from '../../image/EROPEN700.png'
import OneItem from './OneItem';
import OneItem2 from './OneItem2';
import OneItem3 from './OneItem3';
import OneItem4 from './OneItem4';
import OneItem5 from './OneItem5';
import OneItem6 from './OneItem6';
import OneItem7 from './OneItem7';
import Nav from '../nav/Nav';
import { useTranslation } from 'react-i18next';
import i18next from 'i18next';

const Main = () => {
  const { t } = useTranslation()
  const [batafsil, setBatafsil] = useState(true)
  const [ show, setShow] = useState(false)
  useEffect(() => {
    setShow(true)
  }, [])
  const changeLan = (lang) => {
    i18next.changeLanguage(lang)
  }
  useEffect(() => {
  setTimeout(() => {
    const body = document.querySelector("body");
    body.classList.remove("--lock");
    const e = document.querySelector(".about");
    e.classList.contains("--bgAnimate") || e.classList.add("--bgAnimate")
  }, 9500)
},[])



  return (
    <div className='main1' id="mains">
        <Nav changeLan = {changeLan}/>
        <div className='main'>
          <div className={ show? "about --show-item --show --bgAnimate" : "about --show-item --bgAnimate"}>
            <div className="container">
              <div className="about__row">
                <div className= { show ? "about__3d --show" : "about__3d"}>
                  <model-viewer
                    src="./image/scene.gltf"
                    alt="pertin 3d model"
                    auto-rotate
                    disable-zoom
                    camera-controls
                    ios-src="./image/scene.gltf"
                    class="model">
                  </model-viewer>
                  
                </div>
                <div className="about__content --show-item">
                  <h2 className={ show? "about__title --show-item --show" : "about__title --show-item"}>{t("Home.title")}</h2>
                  <a href='#' id='bat' className={show ? " about__link popup-link --show-item --show":"about__link popup-link --show-item"} onClick={() => setBatafsil(!batafsil)}>{t("Home.more")}</a>
                </div>
              </div>
            </div>
          </div>
          <div id="popup" className={batafsil ? "popup" : "popup --open"}>
            <div className="popup__body">
              <div className="popup__content">
                <a href="#" className="popup__close close-popup" onClick={() => setBatafsil(!batafsil)}>X</a>
                <div className="popup__row">
                  <div className="popup__image">
                    <img src={Popup} alt="zaiflik, erkak, man, jinsiy zaiflik, 18+, sex, dori apteka" className='Popup-image' />
                    <div className="popup__text">
                      <h2 className="popup__title">
                        <div className="brand__logo">
                          <img src={Eropen700} alt="EROPEN700" />
                        </div>
                      </h2>
                      <div className="popup__description block --init" data-spollers="" data-one-spoller="">
                        <OneItem />
                        <OneItem2 />
                        <OneItem3 />
                        <OneItem4 />
                        <OneItem5 />
                        <OneItem6 />
                        <OneItem7 />
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
    </div>
  );
};

export default Main;