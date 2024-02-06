import React, {useState, useEffect} from 'react';
import './social.css'
import {FaFacebookF} from "react-icons/fa"
import {FaTelegramPlane} from "react-icons/fa"
import {FaInstagram} from "react-icons/fa"
import { useTranslation } from 'react-i18next';
const Social = () => {
  const { t } = useTranslation()
  const [ show, setShow] = useState(false)
  useEffect(() => {
    setShow(true)
  }, [])
    return (
        <>
            <div className={show ? "main__social social --show-item --show" : "main__social social --show-item"}>
          <div className="social__row">
            <ul className="social__list">
              <li className="social__item">
                <a className="social__link" href="https://www.facebook.com/Eropenuz-101146958816278" target="_blank" rel="noreferrer">
                  <FaFacebookF/></a>
              </li>
              <li className="social__item">
                <a className="social__link" href="https://instagram.com/eropen.uz" target="_blank" rel="noreferrer">
                  <FaInstagram/></a>
              </li>
              <li className="social__item">
                <a className="social__link" href="https://t.me/LbdWoman_Uz" target="_blank" rel="noreferrer"><FaTelegramPlane/></a>
              </li>
            </ul>
          </div>
        </div>
        <div className={ show ? "main__to-order to-order --show-item --show" : "main__to-order to-order --show-item"}>
          <a className="to-order__link" href="#contact">{t('social.send')}</a>
        </div>
        </>
    );
};

export default Social;