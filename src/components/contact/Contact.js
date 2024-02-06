import React, { useState } from 'react';
import Youtube from '../../image/youtube.png'
import Butterfly from '../../image/8.png'
import Butterfly9 from '../../image/9.png'
import Eropen from '../../image/eropen.png'
import './contact.css'

import { useTranslation } from 'react-i18next';


const Contact = () => {
    const [ name, setName ] = useState('')
    const [ phone, setPhone ] = useState('+')
    const [ messages, setMessages ] = useState('')
    const [right, setRight] = useState(false)
    // const [rights, setRights] = useState(false)
    const [right1, setRight1] = useState(false)
    const [right2, setRight2] = useState(false)
    const [success, setSuccess] = useState(false)
    const { t } = useTranslation()
    const sendMessage = (e) => {
        e.preventDefault();
        let utmSource = new URLSearchParams(document.location.search).get('utm_source');
        let utmCampaign = new URLSearchParams(document.location.search).get('utm_campaign');
        let token = "5452032731:AAHr-pD3sRSb-2PkaoY1eoQ5S89cOs11cJI"
        let chatID = "-1001449424479";
        let message = `<b>📬 Eropen.uz</b>%0A<b>👤 Исми: </b><i>${name || '-'}</i>%0A<b>📞 Тел рақами:</b><i>${phone || '-'}</i>%0A<b>📍 Источник:</b> ${utmSource || '-'}%0A<b>📍 Компания:</b> ${utmCampaign || '-'}`;
        let url = `https://api.telegram.org/bot${token}/sendMessage?chat_id=${chatID}&text=${message}&parse_mode=html`;
        let apibot = new XMLHttpRequest();
        apibot.open("GET", url, true);
        setPhone('')
        setName('')
        setMessages('')
        if(name.length > 3 &&phone.length > 12 ) {
            apibot.send();
            setSuccess(true)
            setRight(true)
            setRight1(true)
        }else{
            setRight(false)
            setRight1(false)
        }
    }
    
    const handleName = (e) => {
        setName(e.target.value)
        if(name.length >= 3) {
            setRight(true)
        }
        else {
            setRight(false)
            // setRights(false)
            
        }
      
    }
    const handlePhone = (e) => {
        if(e.target.value.length <= 12) {
            setPhone(e.target.value)
            if(phone.length >= 11) {
                setRight1(true)
            }else {
                setRight1(false)
            }
        }
    }

    const handleMessage = (e) => {
        setMessages(e.target.value)
        if (messages.length >= 4) {
            setRight2(true)
        } else {
            setRight2(false)
        }
    }
   

    return (
        <div className="contact --anim-items" id="contact">
            <div className="container1">
                <div className="contact__form form">
                    <div className="contact__title">
                        <h2 className="form__title">{t('connect.contact')}</h2>
                    </div>
                    <form className="form__body" id="form"  onSubmit={sendMessage}>
                        <div className="form__row">
                            <div className="form__item-row">
                                <div className="form__item">
                                    <input 
                                    className={right ? "form__input --required contact--valid" : "form__input --required"} 
                                    type="text" 
                                    name="firstName" 
                                    id="firstName" 
                                    required value={name} 
                                    onChange={handleName}
                                     />
                                    <label className="form__label" htmlFor="firstName">{t('connect.input1')}</label>
                                    <span className="" id="contact__error-1"></span>
                                </div>
                                <div className="form__item">
                                    <input  
                                        className={right1 ? "form__input --required contact--valid" :  "form__input --required"} 
                                        type="number" 
                                        name="phoneNumber" 
                                        id="phoneNumber" 
                                        required 
                                        maxLength={14}
                                        value={phone}  
                                        onChange={handlePhone} 
                                    />
                                    <label className="form__label" htmlFor="phoneNumber">{t('connect.input2')}
                                    </label>
                                    <span className="" id="contact__error-2"></span>
                                </div>
                            </div>
                            <div className="form__fly">
                                <a href="https://www.youtube.com/watch?v=2pK_EAUdEOY" target="_blank" rel='noreferrer' savefrom_lm_index="2" savefrom_lm="1">
                                    <img src={Youtube} alt="Butterfly" />
                                    <div className="form__you-tube">{t('connect.more')}</div>
                                </a>
                                <span className='savefrom'>
                                <a href="https://en.savefrom.net/161/#url=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3D2pK_EAUdEOY&utm_source=userjs-chrome&utm_medium=extensions&utm_campaign=link_modifier" className='savefrom-icon' rel='noreferrer' target="_blank"></a></span>
                            </div>
                        </div>
                        <div className="form__item">
                            <textarea 
                            className={right2 ? "form__input --required contact--valid" : "form__input --required"} 
                            name="message" 
                            id="message" 
                            cols="30" 
                            rows="3" 
                            value={messages}
                            required 
                            onChange={handleMessage}
                            >
                            </textarea>
                            <label className="form__label" htmlFor="message">{t('connect.textarea')}</label>
                            <span id="contact__error-3"></span>
                        </div>
                        <button href="#contact__popup" className="form__button" type="submit">
                        {t('connect.send')}
                        </button>
                    </form>
                </div>
                <div id="contact__popup" className={success ? "popup --open" : "popup"}>
                    <div className="popup__body">
                        <div className="popup__content contact__content">
                            <a href="/" className="popup__close close-popup" onClick={() => setSuccess(!success)}>X</a>
                            <div className="popup__row">
                                <div className="contact__info">
                                    <p className="contact__greater">
                                        Assalamu Alaykum,
                                        <span className="contact__user-name">{name}</span>
                                    </p>
                                    Murojaatingiz uchun tashakkur! Biz siz bilan 30 soniya
                                    ichida <span className="contact__user-number">{phone}</span> raqamiga
                                    bog'lanamiz!
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div id="contact__popup-error" className="popup">
                    <div className="popup__body">
                        <div className="popup__content contact__content">
                            <a href="https://eropen.uz/#" className="popup__close close-popup">X</a>
                            <div className="popup__row">
                                <p className="contact__info">Xatolik yuz berdi!</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <div className="second">
                <div className="container">
                    <div className="second__row --anim-items">
                        <div className="second__item">
                            <div className="second__text">{t('footer.footerTitle')}</div>
                            <div className="second__text">{t('footer.footerTitle1')}</div>
                            <div className="second__flym">
                                <img src={Butterfly} alt="Butterfly" />
                            </div>
                        </div>
                        <div className="second__item">
                            <a className="second__image" href="https://lbdwoman.uz/" rel='noreferrer' target="_blank">
                                <img src={Eropen} alt="LBD Woman" />
                            </a>
                            <div className="second__flyw">
                                <img src={Butterfly9} alt="Butterfly" />
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <footer className="footer">
                <div className="container2">
                    <div className="footer__row">
                        <p className="footer__text">{t('footer.footerText')}</p>
                    </div>
                </div>
            </footer>
        </div>
    );
};

export default Contact;